import os, re


class Case:
    def __init__(self, path):
        self.path = path
        self.name = os.path.basename(path)
        self.path_nam = self.find_file_by_name(".nam", [])
        self.path_log = self.find_file_by_name(".log", ["runscript"])
        self.path_out = self.find_file_by_name(".out", ["gmon", "mcs.", "mesh."])
        self.powerICP = 0
        self.powerCCP1 = 0
        self.powerCCP2 = 0
        self.dc_bias = eval(self.find_out_parameter("DC BIAS"))
        self.ne_ave = eval(self.find_out_parameter("AVERAGE ELECTRON DENSITY"))
        self.restart = int(self.find_nam_parameter("IRESTART"))
        print(self.restart)


    def find_file_by_name(self, str_match, str_exclude=[]):
        """
        searches case folder for file, matching base on filename

        Args:
            str_match: string to match
            str_exclude: string to exclude

        Returns:
            Path to file or "None" if not found

        """
        # loop over files in case directory
        for name_file in os.listdir(self.path):
            if str_match in name_file.lower() and not any(x in name_file.lower() for x in str_exclude):
                return os.path.join(self.path, name_file)
        return None

    def find_nam_parameter(self, str_match, multiple_values=False):

        with open(self.path_nam, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if re.search(rf' *{str_match} *=', line):
                    if "!" in line:
                        parameter = re.findall(r'= *(.*),\s*\!', line)[0]
                    else:
                        parameter = re.findall(r'= *(.*),\s*$', line)[0]
                    return parameter
        return None

    def find_out_parameter(self, str_match):
        with open(self.path_out, 'r') as f:
            lines = f.readlines()
            for line in lines[::-1]:
                if re.search(rf' *{str_match} *', line):

                    parameter = re.findall(rf' *{str_match} *(.*) ', line)[0]
                    return parameter
        return None


# initialize case list
cases = []

# root directory
dir_root = os.path.abspath("C:\\Users\\flori\\Desktop\\temp2\\ConstVoltage_200_1000")
# dir_root = os.path.abspath("D:\\UIGEL5_D_Florian\Voltage_Waveform_Tailoring\HPEM\ArCF4O2\DarkSpaceGeometry\ConstVoltage_200_1000")
dirs_in_root = os.listdir(dir_root)
dir_source = ""

# find source and case folders
# --------------------------------------------------------------
print("scanning directory...")
for dir_in_root in dirs_in_root:

    dir_sub = os.path.join(dir_root, dir_in_root)
    # check if file:
    string_dir_type = "{0:<17}".format("not recognised")
    if os.path.isfile(dir_sub):
        continue

    # check if directory
    if os.path.isdir(dir_sub):

        # if directory, find out if source or case directory or neither
        # -------------------------------------------------------------------
        paths_in_sub = os.listdir(dir_sub)
        # loop over files and folders
        for path_in_sub in paths_in_sub:

            # check if its a file

            if os.path.isfile(os.path.join(dir_sub, path_in_sub)):
                # check for specific files to determine if it is a source file
                # ---------------------------------------------------------------------
                # if it's a source file label the directory as source directory

                if "icp.f" in path_in_sub.lower() or "icp.for" in path_in_sub.lower():
                    dir_source = os.path.join(path_in_sub)
                    break

                if ".log" in path_in_sub.lower():
                    path_case = os.path.join(dir_sub)
                    cases.append(Case(path_case))
                    break

            # if it is a directory, search directory
            else:
                paths_in_sub_sub = os.listdir(os.path.join( dir_sub, path_in_sub))
                for path_in_sub_sub in paths_in_sub_sub:

                    if ".log" in path_in_sub_sub.lower():
                        path_case = os.path.join(dir_sub, path_in_sub)
                        cases.append(Case(path_case))
                        break

print(f"{len(cases)} case directories found:")
for case in cases:
    print("   " + case.path)


