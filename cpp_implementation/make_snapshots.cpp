#include "wave_analysis.h"     // ensure WaveAnalysis is complete for unique_ptr in headers
#include "fitzhugh_nagumo.h"
#include "finite_difference.h"
#include "wave_loader.h"

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_context.h>

#include <Eigen/Dense>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

// ---------- CSV helpers ----------
static void write_csv_row(std::ofstream& ofs, const Vector& v) {
  for (int i = 0; i < v.size(); ++i) { if (i) ofs << ','; ofs << v[i]; }
  ofs << '\n';
}
static void write_csv_row(std::ofstream& ofs, const std::vector<double>& v) {
  for (size_t i = 0; i < v.size(); ++i) { if (i) ofs << ','; ofs << v[i]; }
  ofs << '\n';
}

// ---------- CVODE RHS (reaction + diffusion on u1) ----------
struct RHSData {
  Eigen::SparseMatrix<double> lap; // periodic Laplacian
  Vector params;                   // FHN parameters
  int N;                           // grid points
};

static int rhs(double /*t*/, N_Vector y, N_Vector ydot, void* user_data) {
  auto* data = static_cast<RHSData*>(user_data);
  const int N = data->N;

  double* y_ptr = NV_DATA_S(y);
  double* d_ptr = NV_DATA_S(ydot);

  // state layout: [u1(0..N-1) | u2(0..N-1)]
  Eigen::Map<Vector> u1(y_ptr, N);
  Eigen::Map<Vector> u2(y_ptr + N, N);
  Eigen::Map<Vector> du1(d_ptr, N);
  Eigen::Map<Vector> du2(d_ptr + N, N);

  // local reaction
  for (int i = 0; i < N; ++i) {
    double xx[2] = {u1[i], u2[i]}, dx[2] = {0.0, 0.0};
    FitzHughNagumoSolver::fhn_local(dx, xx, data->params);
    du1[i] = dx[0];
    du2[i] = dx[1];
  }

  // diffusion on u1 only
  du1 += data->lap * u1;
  return 0;
}

// ---------- CLI ----------
struct CLI {
  std::string fastU, fastP, slowU, slowP, outdir = "./snap_out";
  int gg = 1, n = 200;
  double theta = 0.0, xs = 20.0, A = 0.0, t0 = 0.0, tf = 400.0;
};

static void usage() {
  std::cerr
    << "Usage: make_snapshots --fastU <path> --fastP <path> --slowU <path> --slowP <path>\n"
    << "                      [--gg INT] [--theta VAL] [--xs VAL] [--A VAL]\n"
    << "                      [--n INT] [--t0 VAL] [--tf VAL] [--out DIR]\n";
}

static CLI parse_cli(int argc, char** argv) {
  CLI c;
  for (int i = 1; i < argc; ++i) {
    std::string k = argv[i];
    auto need = [&](const char* name) -> std::string {
      if (i + 1 >= argc) { std::cerr << "Missing value for " << name << "\n"; std::exit(2); }
      return std::string(argv[++i]);
    };
    if (k == "--fastU") c.fastU = need("--fastU");
    else if (k == "--fastP") c.fastP = need("--fastP");
    else if (k == "--slowU") c.slowU = need("--slowU");
    else if (k == "--slowP") c.slowP = need("--slowP");
    else if (k == "--gg") c.gg = std::stoi(need("--gg"));
    else if (k == "--theta") c.theta = std::stod(need("--theta"));
    else if (k == "--xs" || k == "--XS") c.xs = std::stod(need("--xs"));
    else if (k == "--A") c.A = std::stod(need("--A"));
    else if (k == "--n") c.n = std::stoi(need("--n"));
    else if (k == "--t0") c.t0 = std::stod(need("--t0"));
    else if (k == "--tf") c.tf = std::stod(need("--tf"));
    else if (k == "--out") c.outdir = need("--out");
    else { std::cerr << "Unknown arg: " << k << "\n"; usage(); std::exit(2); }
  }
  if (c.fastU.empty() || c.fastP.empty() || c.slowU.empty() || c.slowP.empty()) {
    usage(); std::exit(2);
  }
  if (c.n < 2) { std::cerr << "n must be >= 2\n"; std::exit(2); }
  return c;
}

int main(int argc, char** argv) {
  CLI cli = parse_cli(argc, argv);
  fs::create_directories(cli.outdir);

  // 0) SUNDIALS 6.x context
  SUNContext sunctx;
  if (SUNContext_Create(nullptr, &sunctx) != 0) {
    std::cerr << "SUNContext_Create failed\n"; return 1;
  }

  // 1) solver & states from your wave files
  FitzHughNagumoSolver solver;
  States states = WaveLoader::create_states_from_files(
      cli.fastU, cli.fastP, cli.slowU, cli.slowP, solver);

  // If gg maps to a parameter entry, set it here (adjust index/name accordingly)
  // states.params[PARAM_INDEX_GG] = static_cast<double>(cli.gg);

  // 2) initial condition: u0 + X([xs, theta, A])  (u1 gets +A in the window, u2 unchanged)
  std::vector<double> q = {cli.xs, cli.theta, cli.A};
  Matrix perturb = solver.X(q);     // 2 x N
  Matrix u = states.u0 + perturb;   // 2 x N
  const int N = static_cast<int>(u.cols());

  // 3) pack into CVODE state vector (u1||u2)
  N_Vector y = N_VNew_Serial(2 * N, sunctx);
  if (!y) { std::cerr << "N_VNew_Serial failed\n"; SUNContext_Free(&sunctx); return 1; }
  double* yp = NV_DATA_S(y);
  for (int i = 0; i < N; ++i) yp[i]     = u(0, i);
  for (int i = 0; i < N; ++i) yp[N + i] = u(1, i);

  // 4) set up CVODE — stiff config like your repo
  void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (!cvode_mem) { std::cerr << "CVodeCreate failed\n"; N_VDestroy(y); SUNContext_Free(&sunctx); return 1; }

  RHSData data{ solver.get_fd_operator().get_laplacian(), states.params, N };
  const double t0 = cli.t0, tf = cli.tf;

  if (CVodeInit(cvode_mem, rhs, t0, y) != CV_SUCCESS) {
    std::cerr << "CVodeInit failed\n";
    CVodeFree(&cvode_mem); N_VDestroy(y); SUNContext_Free(&sunctx); return 1;
  }

  // tolerances (use scalar tolerances like FHN solver)
  const double reltol = 1e-6;
  const double abstol = 1e-8;
  CVodeSStolerances(cvode_mem, reltol, abstol);

  CVodeSetUserData(cvode_mem, &data);
  CVodeSetMaxNumSteps(cvode_mem, 100000);      // match FHN solver setting
  CVodeSetMaxStep(cvode_mem, 1.0);             // critical: set max step size like FHN solver
  CVodeSetMaxOrd(cvode_mem, 5);                // BDF max order
  CVodeSetStopTime(cvode_mem, tf);

  // Krylov linear solver
  SUNLinearSolver LS = SUNLinSol_SPGMR(y, PREC_NONE, 0, sunctx);
  if (CVodeSetLinearSolver(cvode_mem, LS, nullptr) != CV_SUCCESS) {
    std::cerr << "CVodeSetLinearSolver failed\n";
    SUNLinSolFree(LS); CVodeFree(&cvode_mem); N_VDestroy(y); SUNContext_Free(&sunctx); return 1;
  }

  // 5) snapshot cadence
  const int K = cli.n;
  const double shot_dt = (tf - t0) / (K - 1);
  double next_t = t0;
  int k = 0;

  // 6) metadata & spatial grid
  {
    std::ofstream meta(cli.outdir + "/meta.json");
    meta << std::fixed << std::setprecision(10)
         << "{\n"
         << "  \"gg\": "    << cli.gg    << ",\n"
         << "  \"theta\": " << cli.theta << ",\n"
         << "  \"xs\": "    << cli.xs    << ",\n"
         << "  \"A\": "     << cli.A     << ",\n"
         << "  \"t0\": "    << t0        << ",\n"
         << "  \"tf\": "    << tf        << ",\n"
         << "  \"n\": "     << K         << ",\n"
         << "  \"N\": "     << N         << "\n"
         << "}\n";
  }
  {
    const auto& grid = solver.get_fd_operator().get_grid(); // Eigen::VectorXd
    std::vector<double> x(grid.data(), grid.data() + grid.size());
    std::ofstream xcsv(cli.outdir + "/x.csv");
    write_csv_row(xcsv, x);
  }

  auto dump = [&](int idx){
    char name[256]; std::snprintf(name, sizeof(name), "%s/snap_%04d.csv", cli.outdir.c_str(), idx);
    std::ofstream s(name);
    Vector u1(N);
    for (int i = 0; i < N; ++i) u1[i] = NV_DATA_S(y)[i]; // first N entries are u1
    write_csv_row(s, u1);
  };

  // 7) initial snapshot
  dump(k++);
  next_t += shot_dt;

  // 8) integrate → dump evenly-spaced snapshots
  double tcur = t0;
  while (tcur < tf - 1e-14 && k < K) {
    double tout = next_t;
    int flag = CVode(cvode_mem, tout, y, &tcur, CV_NORMAL);
    if (flag < 0) {
      std::cerr << "CVode failure at t=" << tcur << " (flag=" << flag << ")\n";
      break;
    }
    dump(k++);
    next_t += shot_dt;
  }

  // 9) cleanup
  SUNLinSolFree(LS);
  CVodeFree(&cvode_mem);
  N_VDestroy(y);
  SUNContext_Free(&sunctx);

  std::cerr << "Wrote " << k << " snapshots to " << cli.outdir << "\n";
  return 0;
}
