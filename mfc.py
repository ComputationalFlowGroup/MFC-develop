#!/usr/bin/env python3

import os
import re
import sys
import shutil
import tarfile
import argparse
import urllib.request


MFC_USE_SUBDIR    = ".mfc"
MFC_CONF_FILEPATH = "mfc.conf.yaml"
MFC_LOCK_FILEPATH = f"{MFC_USE_SUBDIR}/mfc.lock.yaml"


class MFCException(Exception):
    pass


def execute_shell_command_safe(command: str, no_exception: bool = False):
    status = os.system(command)

    if status != 0 and not(no_exception):
        raise MFCException(f'Failed to execute command "{command}".')
    
    return status


def import_module_safe(import_name: str, pip_name: str):
    try:
        globals()[import_name] = __import__(import_name)

        return

    except ImportError:
        pass
    
    prompt = f"The Python package {pip_name} needs to be installed. Would you like to install it now? (from ./python_packages/{pip_name}) (yY/nN): "

    if input(prompt).upper().strip() != "Y":
        raise Exception(f'You requested not to download the Python package {pip_name}. Aborting...')

    if 0 != execute_shell_command_safe(f'cd "python_packages/{pip_name}" && python3 setup.py install --user', no_exception=True):
        raise Exception(f'Failed to automatically install the Python package {pip_name}. Please install it manually using Pip. Aborting...')

    print("\n")

    # Restart mfc.py
    os.execv(sys.argv[0], sys.argv)


def clear_line():
    sys.stdout.write("\033[K")


def file_load_yaml(filepath: str):
    try:
        with open(filepath, "r") as f:
            return yaml.safe_load(f)
    
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f'Failed to load YAML from "{filepath}": {exc}')


def file_dump_yaml(filepath: str, data):
    try:
        with open(filepath, "w") as f:
            yaml.dump(data, f)
    
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f'Failed to dump YAML to "{filepath}": {exc}.')


# TODO: Better solution
def uncompress_archive_to(archive_filepath: str, destination: str):
    archive = tarfile.open(archive_filepath, "r")
    archive.extractall("/".join(destination.rstrip("/").split("/")[:-1]))
    archive.close()

    src = "/".join(destination.rstrip("/").split("/")[:-1])+"/"+archive_filepath.split("/")[-1].replace(".tar.gz", "")

    os.rename(src, destination)


def delete_file_safe(filepath: str):
    if os.path.exists(filepath):
        os.remove(filepath)


def create_file_safe(filepath: str):
    if not os.path.exists(filepath):
        open(filepath, "w").close()


def delete_directory_recursive_safe(directory_path: str):
    if os.path.isdir(directory_path):
        shutil.rmtree(directory_path)


def create_directory_safe(directory_path: str):
    os.makedirs(directory_path, exist_ok=True)


def center_ansi_escaped_text(message: str):
    nCols = shutil.get_terminal_size((80, 20)).columns

    to_escape = [re.escape(colorama.Style.RESET_ALL)]
    for key in dir(colorama.Fore):
        if not callable(getattr(colorama.Fore, key)) and not key.startswith("__"):
            to_escape.append(re.escape(getattr(colorama.Fore, key)))

    longest_string_len = max([len(re.compile("|".join(to_escape), flags=re.DOTALL).sub("", line)) for line in message.splitlines()])

    padding = " "*((nCols - longest_string_len) // 2)

    return "\n".join([f'{padding}{line}{padding}' for line in message.splitlines()])


def clear_print(message, end='\n'):
    clear_line()
    print(message, end=end)


def update_symlink(at: str, to: str):
    if os.path.islink(at):
        os.remove(at)
            
    os.symlink(to, at)


class MFCConf:
    def __init__(self):
        self.data = file_load_yaml(MFC_CONF_FILEPATH)

    def __getitem__(self, key: str, default=None):
        if key not in self.data:
            if default==None:
                raise MFCException(f'MFCConf: Key "{key}" doesn\'t exist.')
            else:
                return default

        return self.data[key]

    def get_target_matches(self, name: str):
        return list(filter(lambda x: x["name"] == name, self["targets"]))

    def does_target_exist(self, name: str):
        return len(self.get_target_matches(name)) > 0
    
    def does_unique_target_exist(self, name: str):
        return len(self.get_target_matches(name)) == 1

    def get_target(self, name: str):
        matches = self.get_target_matches(name)
        
        if len(matches) == 0:
            raise MFCException(f'Failed to retrieve dependency "{name}" in "{obj}" with restrict_cc="{restrict_cc}".')

        if len(matches) > 1:
            raise MFCException(f'More than one dependency to choose from for "{name}", defined in "{obj}" with restrict_cc="{restrict_cc}".')

        return matches[0]

    def get_dependency_names(self, name: str, recursive=False, visited: list = None):
        result: list = []

        if visited == None:
            visited = []

        if name not in visited:
            visited.append(name)

            desc = self.get_target(name)

            for dependency_name in desc.get("depends", []):
                result.append(dependency_name)

                if recursive:
                    result  += self.get_dependency_names(dependency_name, recursive=recursive, visited=visited)
                    visited += result

        return list(set(result))


class MFCLock:
    def __init__(self):
        if not os.path.exists(MFC_LOCK_FILEPATH):
            with open(MFC_LOCK_FILEPATH, 'w') as f:
                f.write("targets: []")

        self.data = file_load_yaml(MFC_LOCK_FILEPATH)

    def __getitem__(self, key: str, default=None):
        if key not in self.data:
            if default==None:
                raise MFCException(f'MFCLock: Key "{key}" doesn\'t exist.')
            else:
                return default

        return self.data[key]

    def get_target_matches(self, name: str, restrict_cc: str = None):
        def peek_filter(e: dict):
            if e["name"] != name:
                return False
            
            if restrict_cc is None:
                return True
            
            return e.get("compiler_configuration", None) == restrict_cc
        
        return list(filter(peek_filter, self["targets"]))

    def does_target_exist(self, name: str, restrict_cc: str = None):
        return len(self.get_target_matches(name, restrict_cc)) > 0
    
    def does_unique_target_exist(self, name: str, restrict_cc: str = None):
        return len(self.get_target_matches(name, restrict_cc)) == 1

    def get_target(self, name: str, restrict_cc: str = None):
        matches = self.get_target_matches(name, restrict_cc)
        
        if len(matches) == 0:
            raise MFCException(f'Failed to retrieve dependency "{name}" in "{obj}" with restrict_cc="{restrict_cc}".')

        if len(matches) > 1:
            raise MFCException(f'More than one dependency to choose from for "{name}", defined in "{obj}" with restrict_cc="{restrict_cc}".')

        return matches[0]

    def save(self):
        file_dump_yaml(MFC_LOCK_FILEPATH, self.data)


class MFCArgs:
    def __init__(self, conf: MFCConf):
        parser = argparse.ArgumentParser(description="Wecome to the MFC master script.", )

        compiler_configuration_names = [e["name"] for e in conf["compilers"]["configurations"]]

        parser.add_argument("--build", action="store_true", help="Build targets.")
        parser.add_argument("--test",  action="store_true", help="Test targets.")
        parser.add_argument("--clean", action="store_true", help="Clean the targets.")
        parser.add_argument("-sc", "--set-current", type=str, choices=compiler_configuration_names, 
                            help="Select a compiler configuration to use when running MFC.")

        compiler_target_names = [e["name"] for e in conf["targets"]]
        parser.add_argument("-t", "--targets", nargs="+", type=str,
                            choices=compiler_target_names, default=["MFC"],
                            help="The space-separated targets you wish to have built.")

        parser.add_argument("-cc", "--compiler-configuration", type=str,
                            choices=compiler_configuration_names, default="release")

        parser.add_argument("-j", "--jobs", metavar="N", type=int,
                            help="Allows for N concurrent jobs.", default=1)

        parser.add_argument("-s", "--scratch", action="store_true",
                            help="Build all targets from scratch.")

        self.data = vars(parser.parse_args())

    def __getitem__(self, key: str, default=None):
        if key not in self.data:
            if default==None:
                raise MFCException(f'MFCArgs: Key "{key}" doesn\'t exist.')
            else:
                return default

        return self.data[key]


class MFC:
    def get_base_path(self, cc: str = None):
        if cc is None:
            cc = self.args["compiler_configuration"]
        
        return f'{self.ROOT_PATH}/{MFC_USE_SUBDIR}/{cc}'

    def get_source_path(self, name: str):
        return f'{self.get_base_path()}/src/{name}'

    def get_build_path(self):
        return f'{self.get_base_path()}/build'

    def get_log_filepath(self, name: str):
        return f'{self.get_base_path()}/log/{name}.log'

    def get_temp_path(self, name: str):
        return f'{self.get_base_path()}/temp/{name}'

    def setup_directories(self):
        create_directory_safe(MFC_USE_SUBDIR)

        for d in ["src", "build", "log", "temp"]:
            for cc in [ cc["name"] for cc in self.conf["compilers"]["configurations"] ]:
                create_directory_safe(f"{MFC_USE_SUBDIR}/{cc}/{d}")
                if d == "build":
                    for build_subdir in ["bin", "include", "lib", "share"]:
                        create_directory_safe(f"{MFC_USE_SUBDIR}/{cc}/{d}/{build_subdir}")        

    def check_environment(self):
        print("|--> Checking for the presence of required command-line utilities...", end='\r')

        required = ["python3", "python3-config", "make", "git"]
        required += self.conf["compilers"]["regular"].values()
        required += self.conf["compilers"]["mpi"].values()

        for index, utility in enumerate(required):
            clear_print(f"|--> {index+1}/{len(required)} Checking for {utility}...", end='\r')

            if shutil.which(utility) is None:
                raise MFCException(
                    f'Failed to find the command line utility "{utility}". Please install it or make it visible.')

        # TODO: MacOS Checks
        if sys.platform == "darwin": # MacOS
            pass
            #mfc_state.conf["compilers"]["mpi"]["fortran"]
        
        clear_print(f"|--> Build environment: Passing. ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})")

    def string_replace(self, dependency_name: str, string: str, recursive=True):
        dep = self.conf.get_target(dependency_name)

        regulars = self.conf["compilers"]["regular"]
        mpis     = self.conf["compilers"]["mpi"]

        compiler_cfg = None
        for c_cfg in self.conf["compilers"]["configurations"]:
            if c_cfg["name"] == self.args["compiler_configuration"]:
                compiler_cfg = c_cfg
                break

        if compiler_cfg is None:
            raise MFCException(
                f'Failed to locate the compiler configuration "{self.args["compiler_configuration"]}".')

        if dep["type"] in ["clone", "download"]:
            install_path = self.get_build_path()
            source_path  = self.get_source_path(dependency_name)
        elif dep["type"] == "source":
            install_path = "ERR_INSTALL_PATH_IS_UNDEFINED"
            source_path  = dep["source"]["source"]
        else:
            raise MFCException(f'Unknown type "{dep["type"]}".')

        flags = compiler_cfg["flags"]

        replace_list = [
            ("${MFC_ROOT_PATH}",     self.ROOT_PATH),
            ("${CONFIGURE_OPTIONS}", f'--prefix="{install_path}"'),
            ("${SOURCE_PATH}",       source_path),
            ("${INSTALL_PATH}",      install_path),
            ("${INSTALL_PATH}",      install_path),
            ("${MAKE_OPTIONS}",      f'-j {self.args["jobs"]}'),
            ("${COMPILER_FLAGS}",    f'CFLAGS="{flags.get("c", "")}" CPPFLAGS="{flags.get("c++", "")}" FFLAGS="{flags.get("fortran", "")}"'),
            ("${REGULAR_COMPILERS}", f'CC="{regulars.get ("c", "")}" CXX="{regulars.get  ("c++", "")}" FC="{regulars.get ("fortran", "")}"'),
            ("${MPI_COMPILERS}",     f'CC="{mpis.get     ("c", "")}" CXX="{mpis.get      ("c++", "")}" FC="{mpis.get     ("fortran", "")}"')
        ]

        for e in replace_list:
            string = string.replace(*e)

        # Combine different assignments to flag variables (CFLAGS, FFLAGS, ...)
        for FLAG_NAME in [ "CFLAGS", "CPPFLAGS", "FFLAGS" ]:
            FLAG_PATTERN = f' {FLAG_NAME}=".*?"'

            matches = re.findall(FLAG_PATTERN, string)

            if len(matches) <= 1:
                continue

            for i in range(len(matches)):
                matches[i] = re.sub(f'^ {FLAG_NAME}="', ' ', matches[i])
                matches[i] = re.sub(r'"$', ' ', matches[i])

            string = re.sub(FLAG_PATTERN, ' ', string, len(matches) - 1)
            string = re.sub(FLAG_PATTERN, f' {FLAG_NAME}="{" ".join(matches)}"', string)                

        # Fetch 
        if recursive:
            for dep2_info in self.conf["targets"]:
                string = string.replace("${" + dep2_info["name"] + ":", "${")
                string = self.string_replace(dep2_info["name"], string, recursive=False)

        return string

    def is_build_satisfied(self, name: str, ignoreIfSource: bool = False):
        # Check if it hasn't been built before
        if not self.lock.does_target_exist(name, self.args["compiler_configuration"]):
            return False
        
        # Retrive CONF & LOCK descriptors
        conf_desc = self.conf.get_target(name)
        lock_desc = self.lock.get_target(name, self.args["compiler_configuration"])

        # Check if it needs updating (LOCK & CONFIG descriptions don't match)
        if conf_desc["type"] != lock_desc["type"]                or \
           conf_desc["type"] == "source" and not(ignoreIfSource) or \
           conf_desc[conf_desc["type"]] != conf_desc[conf_desc["type"]]:
            return False

        # Check if any of its dependencies needs updating
        for dependency_name in self.conf.get_dependency_names(name, recursive=True):
            if not self.is_build_satisfied(dependency_name, ignoreIfSource=ignoreIfSource):
                return False

        # Check for "scratch" flag
        if self.args["scratch"]:
            return False

        return True

    def build_target__clean_previous(self, name: str):
        if not self.lock.does_unique_target_exist(name, self.args["compiler_configuration"]):
            return
        
        conf_desc = self.conf.get_target(name)
        lock_desc = self.lock.get_target(name, self.args["compiler_configuration"])
        
        if ((    conf_desc["type"] != lock_desc["type"]
             and lock_desc["type"] in ["clone", "download"]
            ) or (self.args["scratch"])):
            delete_directory_recursive_safe(f'{self.ROOT_PATH}/{MFC_USE_SUBDIR}/{lock_desc["compiler_configuration"]}/src/{name}')

    def build_target__fetch(self, name: str, logfile):
        conf = self.conf.get_target(name)
        
        if conf["type"] in ["clone", "download"]:
            if conf["type"] == "clone":
                lock_matches = self.lock.get_target_matches(name, self.args["compiler_configuration"])
                
                if ((   len(lock_matches)    == 1
                    and conf["clone"]["git"] != self.lock.get_target(name, self.args["compiler_configuration"])["clone"]["git"])
                    or (self.args["scratch"])):
                    clear_print(f'|--> Package {name}: GIT repository changed. Updating...', end='\r')

                    delete_directory_recursive_safe(self.get_source_path(name))

                if not os.path.isdir(self.get_source_path(name)):
                    clear_print(f'|--> Package {name}: Cloning repository...', end='\r')
                    
                    execute_shell_command_safe(
                        f'git clone --recursive "{conf["clone"]["git"]}" "{self.get_source_path(name)}" >> "{logfile.name}" 2>&1')

                clear_print(f'|--> Package {name}: Checking out {conf["clone"]["hash"]}...', end='\r')

                execute_shell_command_safe(
                    f'cd "{self.get_source_path(name)}" && git checkout "{conf["clone"]["hash"]}" >> "{logfile.name}" 2>&1')
            elif conf["type"] == "download":
                clear_print(f'|--> Package {name}: Removing previously downloaded version...', end='\r')

                delete_directory_recursive_safe(self.get_source_path(name))

                download_link = conf[conf["type"]]["link"].replace("${VERSION}", conf[conf["type"]]["version"])
                filename = download_link.split("/")[-1]

                clear_print(f'|--> Package {name}: Downloading source...', end='\r')

                create_directory_safe(self.get_temp_path(name))
                
                download_path = f'{self.get_temp_path(name)}/{filename}'                
                urllib.request.urlretrieve(download_link, download_path)

                clear_print(f'|--> Package {name}: Uncompressing archive...', end='\r')

                uncompress_archive_to(download_path,
                                      f'{self.get_source_path(name)}')

                os.remove(download_path)
        elif conf["type"] == "source":
            pass
        else:
            raise MFCException(f'Dependency type "{conf["type"]}" is unsupported.')

    def build_target__build(self, name: str, logfile):
        conf = self.conf.get_target(name)

        if conf["type"] not in ["clone", "download", "source"]:
            raise MFCException(f'Unknown target type "{conf["type"]}".')
        
        for command in conf["build"]:
            command = self.string_replace(name, f"""\
cd "${{SOURCE_PATH}}" && \
PYTHON="python3" PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --includes) $(python3-config --libs)" \
bash -c '{command}' >> "{logfile.name}" 2>&1""")

            logfile.write(f'\n--- ./mfc.py ---\n{command}\n--- ./mfc.py ---\n\n')

            clear_print(f'|--> Package {name}: Building (Logging to {logfile.name})...', end='\r')

            execute_shell_command_safe(command)    

    def build_target__update_lock(self, name: str):
        conf = self.conf.get_target(name)

        clear_print(f'|--> Package {name}: Updating lock file...', end='\r')

        new_entry = dict(conf, **{"compiler_configuration": self.args["compiler_configuration"]})

        if len(self.lock.get_target_matches(name, self.args["compiler_configuration"])) == 0:
            self.lock["targets"].append(new_entry)
        else:
            for index, dep in enumerate(self.lock["targets"]):
                if dep["name"] == name:
                    self.lock["targets"][index] = new_entry

        self.lock.save()

    def build_target(self, name: str):
        # Check if it needs to be (re)built
        if self.is_build_satisfied(name):
            clear_print(f'|--> Package {name}: Nothing to do ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})')
            return False
        
        # Build its dependencies
        for dependency_names in self.conf.get_dependency_names(name, recursive=False):
            self.build_target(dependency_names)

        clear_print(f'|--> Package {name}: Preparing build...', end='\r')

        with open(self.get_log_filepath(name), "w") as logfile:
            self.build_target__clean_previous(name)          # Clean any old build artifacts
            self.build_target__fetch         (name, logfile) # Fetch Source Code
            self.build_target__build         (name, logfile) # Build
            self.build_target__update_lock   (name)          # Update LOCK

        clear_print(f'|--> Package {name}: Done. ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})')

        return True

    def print_header(self):
        print(center_ansi_escaped_text(f"""{colorama.Fore.BLUE}
     ___            ___          ___     
    /__/\          /  /\        /  /\    
   |  |::\        /  /:/_      /  /:/    
   |  |:|:\      /  /:/ /\    /  /:/     
 __|__|:|\:\    /  /:/ /:/   /  /:/  ___ 
/__/::::| \:\  /__/:/ /:/   /__/:/  /  /\ 
\  \:\~~\__\/  \  \:\/:/    \  \:\ /  /:/
 \  \:\         \  \::/      \  \:\  /:/ 
  \  \:\         \  \:\       \  \:\/:/  
   \  \:\         \  \:\       \  \::/   
    \__\/          \__\/        \__\/    
{colorama.Style.RESET_ALL}"""))

        print(center_ansi_escaped_text("""
+----------------------------------+
|    Multi-component Flow Code     |
|                                  |
| https://github.com/MFlowCode/MFC |
+----------------------------------+

"""))

    def test_target(self, name: str):
        if not self.is_build_satisfied(name, ignoreIfSource=True):
            raise MFCException(f"Can't test {name} because its build isn't satisfied.")
        
        with open(self.get_log_filepath(name), "w") as logfile:
            for command in self.conf.get_target(name)["test"]:
                command = self.string_replace(name, f"""\
cd "${{SOURCE_PATH}}" && \
bash -c '{command}' >> "{logfile.name}" 2>&1""")

                logfile.write(f'\n--- ./mfc.py ---\n{command}\n--- ./mfc.py ---\n\n')

                execute_shell_command_safe(command)

                clear_print(f'|--> Package {name}: Testing (Logging to {logfile.name})...', end='\r')

        clear_print(f'|--> Package {name}: Tested. ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})')

    def clean_target(self, name: str):
        delete_directory_recursive_safe(self.get_source_path(name))
        delete_file_safe(self.get_log_filepath(name))

        self.lock["targets"] = list(filter(lambda x: x["name"] != name, self.lock["targets"]))
        self.lock.save()

    def __init__(self):
        self.ROOT_PATH = os.path.dirname(os.path.realpath(__file__))

        self.conf = MFCConf()
        self.setup_directories()
        self.lock = MFCLock()
        self.args = MFCArgs(self.conf)

        self.print_header()
        self.check_environment()

        if self.args["set_current"] is not None:
            update_symlink(f"{MFC_USE_SUBDIR}/___current___", self.get_base_path(self.args["set_current"]))

        # Update symlink to current build
        if self.args["build"]:
            update_symlink(f"{MFC_USE_SUBDIR}/___current___", self.get_base_path())

        for target_name in [ x["name"] for x in self.conf["targets"] ]:
            if target_name in self.args["targets"]:
                if self.args["build"]: self.build_target(target_name)
                if self.args["test"]:  self.test_target(target_name)
                if self.args["clean"]: self.clean_target(target_name)

        self.lock.save()


def main():
    execute_shell_command_safe("git submodule update --init --recursive", no_exception=True)

    for module in [("yaml", "pyyaml"), ("colorama", "colorama")]:
        import_module_safe(*module)
    
    try:
        colorama.init()

        mfc = MFC()
    except (MFCException, Exception) as exc:
        print(f"{colorama.Fore.RED}|--> {str(exc)}{colorama.Style.RESET_ALL}")
        exit(1)


if __name__ == "__main__":
    main()

