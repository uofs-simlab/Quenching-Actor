#!/usr/bin/env python3
"""
Script to analyze and compare job completion data between bracket-enabled and bracket-disabled runs.
Creates a CSV file with job parameters, results, execution times, and differences.
"""

import re
import csv
import sys
from typing import Dict, List, Tuple, Optional

def parse_job_completion_line(line: str) -> Optional[Dict]:
    """Parse a job completion line to extract parameters and results."""
    pattern = r'Job completed in ([\d.]+)s: gg=(\d+), theta=(-?\d+), xs=([\d.]+), n=(\d+), usmin=(-?[\d.]+), usmax=(-?[\d.]+), result=(-?[\d.]+)'
    match = re.search(pattern, line)
    
    if match:
        return {
            'runtime': float(match.group(1)),
            'gg': int(match.group(2)),
            'theta': int(match.group(3)),
            'xs': float(match.group(4)),
            'n': int(match.group(5)),
            'usmin': float(match.group(6)),
            'usmax': float(match.group(7)),
            'result': float(match.group(8))
        }
    return None

def create_job_key(job: Dict) -> Tuple:
    """Create a unique key for a job based on its parameters."""
    return (job['gg'], job['theta'], job['xs'])

def parse_job_file(filename: str) -> Tuple[List[Dict], bool]:
    """Parse a job file and return list of job completions and bracket status."""
    jobs = []
    bracket_enabled = False
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                # Check bracket status
                if 'bracket_enabled=' in line:
                    bracket_enabled = 'bracket_enabled=true' in line
                
                # Parse job completion lines
                job = parse_job_completion_line(line)
                if job:
                    jobs.append(job)
                    
    except FileNotFoundError:
        print(f"Error: File {filename} not found")
        sys.exit(1)
        
    return jobs, bracket_enabled

def main():
    # File paths
    update_file = '/globalhome/tus210/HPC/quenchin_actor/julia_quench_exhustive_update-3717346.out'
    no_update_file = '/globalhome/tus210/HPC/quenchin_actor/julia_quench_exhustive-3720443.out'
    output_file = '/globalhome/tus210/HPC/quenchin_actor/job_comparison_analysis.csv'
    
    print("Parsing job files...")
    
    # Parse both files
    update_jobs, update_bracket_enabled = parse_job_file(update_file)
    no_update_jobs, no_update_bracket_enabled = parse_job_file(no_update_file)
    
    print(f"Update file: {len(update_jobs)} jobs, bracket_enabled={update_bracket_enabled}")
    print(f"No-update file: {len(no_update_jobs)} jobs, bracket_enabled={no_update_bracket_enabled}")
    
    # Create dictionaries with job keys for easy lookup
    update_dict = {create_job_key(job): job for job in update_jobs}
    no_update_dict = {create_job_key(job): job for job in no_update_jobs}
    
    # Find common jobs
    common_keys = set(update_dict.keys()) & set(no_update_dict.keys())
    update_only_keys = set(update_dict.keys()) - set(no_update_dict.keys())
    no_update_only_keys = set(no_update_dict.keys()) - set(update_dict.keys())
    
    print(f"Common jobs: {len(common_keys)}")
    print(f"Update-only jobs: {len(update_only_keys)}")
    print(f"No-update-only jobs: {len(no_update_only_keys)}")
    
    # Prepare CSV data
    csv_data = []
    
    # Add common jobs with comparison
    for key in sorted(common_keys):
        update_job = update_dict[key]
        no_update_job = no_update_dict[key]
        
        runtime_diff = update_job['runtime'] - no_update_job['runtime']
        percent_diff = (runtime_diff / no_update_job['runtime']) * 100 if no_update_job['runtime'] > 0 else 0
        
        csv_data.append({
            'gg': update_job['gg'],
            'theta': update_job['theta'],
            'xs': update_job['xs'],
            'n_update': update_job['n'],
            'n_no_update': no_update_job['n'],
            'usmin_update': update_job['usmin'],
            'usmax_update': update_job['usmax'],
            'usmin_no_update': no_update_job['usmin'],
            'usmax_no_update': no_update_job['usmax'],
            'result_update': update_job['result'],
            'result_no_update': no_update_job['result'],
            'runtime_update': update_job['runtime'],
            'runtime_no_update': no_update_job['runtime'],
            'runtime_diff': runtime_diff,
            'runtime_percent_diff': percent_diff,
            'bracket_enabled_update': update_bracket_enabled,
            'bracket_enabled_no_update': no_update_bracket_enabled,
            'job_status': 'common'
        })
    
    # Add update-only jobs
    for key in sorted(update_only_keys):
        job = update_dict[key]
        csv_data.append({
            'gg': job['gg'],
            'theta': job['theta'],
            'xs': job['xs'],
            'n_update': job['n'],
            'n_no_update': '',
            'usmin_update': job['usmin'],
            'usmax_update': job['usmax'],
            'usmin_no_update': '',
            'usmax_no_update': '',
            'result_update': job['result'],
            'result_no_update': '',
            'runtime_update': job['runtime'],
            'runtime_no_update': '',
            'runtime_diff': '',
            'runtime_percent_diff': '',
            'bracket_enabled_update': update_bracket_enabled,
            'bracket_enabled_no_update': no_update_bracket_enabled,
            'job_status': 'update_only'
        })
    
    # Add no-update-only jobs
    for key in sorted(no_update_only_keys):
        job = no_update_dict[key]
        csv_data.append({
            'gg': job['gg'],
            'theta': job['theta'],
            'xs': job['xs'],
            'n_update': '',
            'n_no_update': job['n'],
            'usmin_update': '',
            'usmax_update': '',
            'usmin_no_update': job['usmin'],
            'usmax_no_update': job['usmax'],
            'result_update': '',
            'result_no_update': job['result'],
            'runtime_update': '',
            'runtime_no_update': job['runtime'],
            'runtime_diff': '',
            'runtime_percent_diff': '',
            'bracket_enabled_update': update_bracket_enabled,
            'bracket_enabled_no_update': no_update_bracket_enabled,
            'job_status': 'no_update_only'
        })
    
    # Write CSV file
    fieldnames = [
        'gg', 'theta', 'xs',
        'n_update', 'n_no_update',
        'usmin_update', 'usmax_update', 'usmin_no_update', 'usmax_no_update',
        'result_update', 'result_no_update',
        'runtime_update', 'runtime_no_update', 'runtime_diff', 'runtime_percent_diff',
        'bracket_enabled_update', 'bracket_enabled_no_update',
        'job_status'
    ]
    
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_data)
    
    print(f"\nCSV file created: {output_file}")
    print(f"Total rows: {len(csv_data)}")
    
    # Print some statistics
    common_jobs = [row for row in csv_data if row['job_status'] == 'common']
    if common_jobs:
        runtime_diffs = [float(row['runtime_diff']) for row in common_jobs if row['runtime_diff'] != '']
        if runtime_diffs:
            avg_diff = sum(runtime_diffs) / len(runtime_diffs)
            print(f"\nRuntime Analysis:")
            print(f"Average runtime difference (update - no_update): {avg_diff:.3f}s")
            print(f"Min runtime difference: {min(runtime_diffs):.3f}s")
            print(f"Max runtime difference: {max(runtime_diffs):.3f}s")
            
            faster_update = sum(1 for d in runtime_diffs if d < 0)
            faster_no_update = sum(1 for d in runtime_diffs if d > 0)
            print(f"Jobs where update was faster: {faster_update}")
            print(f"Jobs where no-update was faster: {faster_no_update}")

if __name__ == "__main__":
    main()
