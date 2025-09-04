#ifndef _CONFIG_H
#define _CONFIG_H

// include the necessary libraries
#include "caf/all.hpp"
#include "caf/io/all.hpp"
#include <chrono>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <queue>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <map>
#include <tuple>
#include <cmath>
#include <cstdio>


using namespace caf;
using namespace std;

class config : public caf::actor_system_config
{
public:
    uint16_t port = 0;
    std::string host = "localhost";
    bool server_mode = false;
    bool enable_bracket = false;
    bool enable_early_termination = false;  
    config()
    {
        opt_group{custom_options_, "global"}
            .add(port, "port,p", "set port")
            .add(host, "host,H", "set host (ignored in server mode)")
            .add(server_mode, "server-mode,s", "enable server mode")
            .add(enable_bracket, "enable-bracket,b", "enable bracket optimization")
            .add(enable_early_termination, "enable-early-termination,e", "enable Sundials early termination");
    }
};

#endif