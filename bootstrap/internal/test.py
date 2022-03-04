#!/usr/bin/env python3

import os
import re
import copy
import hashlib
import colorama
import subprocess
import dataclasses

from pathlib import Path

import internal.common      as common
import internal.treeprint   as treeprint


@dataclasses.dataclass
class Case:
    name:       str
    parameters: list

    def __init__(self, data: dict) -> None:
        self.name       = data.get("name")
        self.parameters = {}
        
        for p in list(data.get("parameters").items()):
            self.parameters[p[0]] = p[1]

    def get_keys(self):
        keys = []
        for param in self.parameters:
            keys.append(param.name)
        
        return keys

    def has_parameter(self, key: str):
        return key in self.get_keys()

    def __getitem__(self, key: str) -> str:
        if key not in self.parameters:
            raise common.MFCException(f"Case: Parameter {key} does not exist.")
        
        return self.parameters[key]

    def __setitem__(self, key: str, val: str):
        self.parameters[key] = val

    def create_case_dict_str(self) -> str: 
        result: str = "{\n"

        for key,val in self.parameters.items():
            result = f'{result}\t"{key}": "{val}",\n'

        return result + "}"


@dataclasses.dataclass
class Test:
    case: Case

    def __init__(self, data: dict) -> None:
        self.case = data.get("case", {})


Tend = 0.25
Nt   = 500
mydt = Tend/(1.*Nt)

BASE_CASE = Case({
    "name": "Base Case",
    "parameters": {
        'case_dir'                     : '\'.\'',
        'run_time_info'                : 'F',
        'nodes'                        : 1,
        'ppn'                          : 1,
        'queue'                        : 'normal',
        'walltime'                     : '24:00:00',
        'mail_list'                    : '',
        'm'                            : 0,       
        'n'                            : 0,        
        'p'                            : 0,        
        'dt'                           : mydt,     
        't_step_start'                 : 0,        
        't_step_stop'                  : int(Nt+1),
        't_step_save'                  : int(Nt),  
        'num_patches'                  : 2,     
        'model_eqns'                   : 2,     
        'alt_soundspeed'               : 'F',   
        'num_fluids'                   : 1,     
        'adv_alphan'                   : 'T',   
        'mpp_lim'                      : 'F',   
        'mixture_err'                  : 'F',   
        'time_stepper'                 : 3,     
        'weno_vars'                    : 2,     
        'weno_order'                   : 5,     
        'weno_eps'                     : 1.E-16,
        'mapped_weno'                  : 'T',
        'null_weights'                 : 'F',
        'mp_weno'                      : 'F',
        'riemann_solver'               : 2,
        'wave_speeds'                  : 1,
        'avg_state'                    : 2,
        'format'                       : 1,
        'precision'                    : 2,
        'prim_vars_wrt'                :'T',
        'parallel_io'                  :'F',

        'patch_icpp(1)%pres'           : 1.0,
        'patch_icpp(1)%alpha_rho(1)'   : 1.E+00,
        'patch_icpp(1)%alpha(1)'       : 1.,
        
        'patch_icpp(2)%pres'           : 0.1,
        'patch_icpp(2)%alpha_rho(1)'   : 0.125E+00,
        'patch_icpp(2)%alpha(1)'       : 1.,
        'fluid_pp(1)%gamma'            : 1.E+00/(1.4-1.E+00),
        'fluid_pp(1)%pi_inf'           : 0.0,
    }
})

class MFCTest:
    def __init__(self, bootstrap):
        self.bootstrap = bootstrap

        # Aliases
        self.tree = self.bootstrap.tree
        self.args = self.bootstrap.args

    def get_test_params(self):
        all_run_params = []

        for dimInfo in [ (["x"],           {'m': 299},                  {"geometry": 1}),
                         (["x", "y"],      {'m': 49, 'n': 39},          {"geometry": 3}),
                         (["x", "y", "z"], {'m': 39, 'n': 29, 'p': 19}, {"geometry": 9}) ]:
            dimParams = {**dimInfo[1]}

            for dimCmp in dimInfo[0]:
                dimParams[f"{dimCmp}_domain%beg"] = 0.E+00
                dimParams[f"{dimCmp}_domain%end"] = 1.E+00

                dimParams[f"bc_{dimCmp}%beg"] = -3
                dimParams[f"bc_{dimCmp}%end"] = -3

            for patchID in [1,2]:
                dimParams[f"patch_icpp({patchID})%geometry"] = dimInfo[2].get("geometry")
                
                if "x" in dimInfo[0]:
                    dimParams[f"patch_icpp({1})%x_centroid"] = 0.25
                    dimParams[f"patch_icpp({2})%x_centroid"] = 0.75
                    dimParams[f"patch_icpp({patchID})%length_x"] = 0.5
                    dimParams[f"patch_icpp({patchID})%vel(1)"] = 0.0
                
                if "y" in dimInfo[0]:
                    dimParams[f"patch_icpp({patchID})%y_centroid"] = 0.5
                    dimParams[f"patch_icpp({patchID})%length_y"]   = 1
                    dimParams[f"patch_icpp({patchID})%vel(2)"]     = 0.0
                
                if "z" in dimInfo[0]:
                    dimParams[f"patch_icpp({patchID})%z_centroid"] = 0.5
                    dimParams[f"patch_icpp({patchID})%length_z"]   = 1
                    dimParams[f"patch_icpp({patchID})%vel(3)"]     = 0.0

            for weno_order in [3, 5]:
                for mapped_weno, mp_weno in [('F', 'F'), ('T', 'F'), ('F', 'T')]:
                    if not (mp_weno == 'T' and weno_order != 5):
                        all_run_params.append({**dimParams, **{'weno_order': weno_order, 'mapped_weno': mapped_weno, 'mp_weno': mp_weno}})

            for riemann_solver in [1, 2]:
                # FIXME: alt_soundspeed not supported for a single fluid
                #all_run_params.append({**dimParams, **{'weno_order': weno_order, 'riemann_solver': riemann_solver, 'alt_soundspeed': 'T'}})
                all_run_params.append({**dimParams, **{'weno_order': weno_order, 'riemann_solver': riemann_solver, 'mixture_err':    'T'}})
                # FIXME: mpp_lim not supported for a single fluid
                #all_run_params.append({**dimParams, **{'weno_order': weno_order, 'riemann_solver': riemann_solver, 'mpp_lim':        'T'}})
                all_run_params.append({**dimParams, **{'weno_order': weno_order, 'riemann_solver': riemann_solver, 'avg_state':      1}})
                all_run_params.append({**dimParams, **{'weno_order': weno_order, 'riemann_solver': riemann_solver, 'wave_speeds':    2}})

                # TODO: num_comp

            all_run_params.append({**dimParams, **{'ppn': 2}})
        
        return all_run_params

    def test(self):
        self.tree.print(f"Testing mfc")
        self.tree.indent()

        if self.args["generate"]:
            common.delete_directory_recursive(common.MFC_TESTDIR)
            common.create_directory(common.MFC_TESTDIR)

        for target in ["pre_process", "simulation"]:
            if not self.bootstrap.is_build_satisfied(target):
                self.tree.print(f"{target} needs (re)building...")
                self.bootstrap.build_target(f"{target}")
        
        all_test_params = self.get_test_params()
        for i, run_params in enumerate(all_test_params):
            self.tree.print_progress(f"Running test #{i+1} - {self.get_case_dir_name(run_params)}", i+1, len(all_test_params))
            self.handle_case(i, run_params)

        common.clear_line()
        self.tree.print(f"Tested. ({colorama.Fore.GREEN}SUCCESS{colorama.Style.RESET_ALL})")
        self.tree.unindent()

    def get_case_dir_name(self, mods: dict):
        return hashlib.sha1("".join([f"{x[0]}{x[1]}" for x in mods.items()]).encode()).hexdigest()[:20]

    def get_case_dir(self, mods: dict):
        return f"{common.MFC_TESTDIR}/{self.get_case_dir_name(mods)}"

    def get_case_from_mods(self, mods: dict):
        case = copy.deepcopy(BASE_CASE)

        for key, val in mods.items():
            case[key] = val
        
        return case

    def create_case_dir(self, mods: dict):
        case     = self.get_case_from_mods(mods)
        case_dir = self.get_case_dir(mods)

        content = f"""\
#!/usr/bin/env python3

from pathlib import Path
from sys     import path

path.insert(0, f"{{Path(__file__).parent.resolve()}}/../../src/master_scripts")

# Let Python find MFC's module
from m_python_proxy import f_execute_mfc_component

case_dict = {case.create_case_dict_str()}

f_execute_mfc_component('pre_process', case_dict, '..', 'serial')
f_execute_mfc_component('simulation',  case_dict, '..', 'serial')

"""

        common.create_directory(case_dir)

        common.file_write(f"{case_dir}/input.py", content)

    def golden_file_compare_match(self, truth: str, candidate: str):
        for candidate_line in candidate.splitlines():
            if candidate_line == "":
                continue

            file_subpath: str = candidate_line.split(' ')[0]
            
            line_trusted: str = ""
            for l in truth.splitlines():
                if l.startswith(file_subpath):
                    line_trusted = l
                    break
            
            if len(line_trusted) == 0:
                continue
            
            numbers_cand  = [ float(x) for x in candidate_line.strip().split(' ')[1:] ]
            numbers_trust = [ float(x) for x in line_trusted.strip().split(' ')[1:]   ]

            # Different amount of spaces, means that there are more entires in one than in the other
            if len(numbers_cand) != len(numbers_trust):
                return (False, "Variable count didn't match.")

            # check values one by one
            for i in range(len(numbers_cand)):
                abs_delta = abs(numbers_cand[i]-numbers_trust[i])
                rel_diff  = abs(abs_delta/numbers_trust[i]) if numbers_trust[i] != 0 else 0
                if    (abs_delta > 1e-12 and rel_diff > 1e-12) \
                   or (numbers_cand[i] * numbers_trust[i] < 0): # Check if sign matches
                    percent_diff = rel_diff*100
                    return (False, f"Error margin is too high for the value #{i+1} in {file_subpath}: ~{round(percent_diff, 5)}% (~{round(abs_delta, 5)}).")

        # Both tests gave the same results within an acceptable tolerance
        return (True, "")

    def get_test_summary(self, mods: dict):
        return "".join([f"{str(x[0]).split('_')[0][:4]}-{str(x[1])[:4]}_" for x in mods.items()])[:-1]

    def handle_case(self, testID, parameters: dict):
        self.create_case_dir(parameters)
        
        def on_test_errror(msg: str = "", term_out: str = ""):
            common.clear_line()
            self.tree.print(f"Test #{testID}: Failed! ({colorama.Fore.RED}FAILURE{colorama.Style.RESET_ALL})")
            if msg != "":
                self.tree.print(msg)
            
            common.file_write(f"{common.MFC_TESTDIR}/failed_test.txt", f"""\
(1/3) Test #{testID}:
  - Summary:  {self.get_test_summary(parameters)}
  - Location: {self.get_case_dir(parameters)}
  - Error:    {msg}

(2/3) Test case:
{self.get_case_from_mods(parameters).create_case_dict_str()}

(3/3) Terminal output:
{term_out}
""")

            self.tree.print(f"Please read {common.MFC_TESTDIR}/failed_test.txt for more information.")
            raise common.MFCException("Testing failed (view above).")

        cmd = subprocess.run(f"cd '{self.get_case_dir(parameters)}' && python3 input.py",
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, shell=True)        
        common.file_write(f"{self.get_case_dir(parameters)}/out.txt", cmd.stdout)
        
        if cmd.returncode != 0:
            on_test_errror("MFC Execution Failed.", cmd.stdout)

        pack = self.pack_case_output(parameters)
        common.file_write(f"{self.get_case_dir(parameters)}/pack.txt", pack)

        golden_filepath = f"{self.get_case_dir(parameters)}/golden.txt"

        if self.args["generate"]:
            common.delete_file(golden_filepath)
            common.file_write(golden_filepath, pack)
        
        if not os.path.isfile(golden_filepath):
            common.clear_line()
            on_test_errror("Golden file doesn't exist! To generate golden files, use the '-g' flag.", cmd.stdout)
        
        golden_file_content = common.file_read(golden_filepath)
        bSuccess, errorMsg  = self.golden_file_compare_match(golden_file_content, pack)
        if not bSuccess:
            on_test_errror(errorMsg, cmd.stdout)

    def pack_case_output(self, params: dict):
        result: str = ""

        case_dir = self.get_case_dir(params)
        D_dir    = f"{case_dir}/D/"

        for filepath in list(Path(D_dir).rglob("*.dat")):
            file_content   = common.file_read(filepath)
            short_filepath = str(filepath).replace(f'{case_dir}/', '')
            
            result += f"{short_filepath} " + re.sub(r' +', ' ', file_content.replace('\n', ' ')).strip() + '\n'
        
        return result
