import os, re
from Movie2XT import movie2xt
import plotting


class Cases:
    def __init__(self):
        self.cases = []

    def __getitem__(self, i):
        return self.cases[i]

    def add_case(self, case):
        self.cases.append(case)

    def __len__(self):
        return len(self.cases)

    def get_value_pair(self, name_var1, name_var2, custom_waveform_only=False):
        """
        extracts value pair from cases and sorts by first value
        Args:
            name_var1 (str): name of first variable (x)
            name_var2 (str): name of first variable (y)
            custom_waveform_only (bool): only consider cases that contain custom waveforms

        Returns: nested list with sorted variables (x and y)

        """
        var1 = []
        var2 = []
        for case in self.cases:
            # check if custom waveform case:
            if custom_waveform_only and not case.contains_custom:
                continue
            var1.append(eval("case."+name_var1))
            var2.append(eval("case."+name_var2))
            # sort lists based on var1
            var1, var2 = (list(t) for t in zip(*sorted(zip(var1, var2))))
        return [var1, var2]


class Case:
    def __init__(self, path):
        self.path = path
        self.name = os.path.basename(path)
        self.path_nam = self.find_file_by_name(".nam", [])
        self.path_log = self.find_file_by_name(".log", ["runscript"])
        self.path_out = self.find_file_by_name(".out", ["gmon", "mcs.", "mesh."])
        self.path_movie1 = self.find_file_by_name("movie1.pdt")
        self.powerICP = 0
        self.powerCCP1 = 0
        self.powerCCP2 = 0
        self.metal_labels = list(eval(self.find_nam_parameter("CMETAL")))
        self.final_voltages = self.get_final_voltages()
        self.rfpnorma = list(eval(self.find_nam_parameter("RFPNORMA")))
        self.dc_bias = eval(self.find_out_parameter("DC BIAS"))
        self.ne_ave = eval(self.find_out_parameter("AVERAGE ELECTRON DENSITY"))
        self.restart = int(self.find_nam_parameter("IRESTART"))
        self.icustom = list(eval(self.find_nam_parameter("ICUSTOM")))
        self.contains_custom = self.icustom != [0] * len(self.icustom)
        self.cwaveform_phase = eval(self.find_nam_parameter("CUSTOM_PHASE"))[1]
        self.rffac = int(float(self.find_nam_parameter("RFFAC")))
        self.freq = float(self.find_nam_parameter("FREQ"))


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
                        parameter = re.findall(r'= *(.*)(?=,*\s*\!)', line)[0]
                    else:
                        parameter = re.findall(r'= *(.*)(?=,*\s*$)', line)[0]
                    # remove leading and trailing whitespaces
                    parameter = parameter.strip()
                    # remove last ',' symbol if present
                    if parameter[-1] == ",":
                        parameter = parameter[:-1]
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

    def get_final_voltages(self):
        voltages = []
        for num_metal, metal in enumerate(self.metal_labels):
            parameter = 0
            with open(self.path_out, 'r') as f:
                lines = f.readlines()
                for line in lines[::-1]:
                    regex = rf' +{metal} +\d\.\d\d\dE\+\d\d  +(\d\.\d\d\dE\+\d\d) +'
                    if re.search(regex, line):
                        parameter = eval(re.findall(regex, line)[0])
                        break
            voltages.append(parameter)
        return voltages

# initialize case list
cases = Cases()

# root directory
# dir_root = os.path.abspath("C:\\Users\\flori\\Desktop\\temp2\\ConstVoltage_200_1000")
#dir_root = os.path.abspath("D:\\UIGEL5_D_Florian\Voltage_Waveform_Tailoring\HPEM\ArCF4O2\DarkSpaceGeometry\ConstVoltage_200_1000")
dir_root = os.path.abspath("D:\\UIGEL5_D_Florian\\Voltage_Waveform_Tailoring\\HPEM\\ArCF4O2\\DarkSpaceGeometry\\2000_2000")
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
                    cases.add_case(Case(path_case))
                    break

            # if it is a directory, search directory
            else:
                paths_in_sub_sub = os.listdir(os.path.join( dir_sub, path_in_sub))
                for path_in_sub_sub in paths_in_sub_sub:

                    if ".log" in path_in_sub_sub.lower():
                        path_case = os.path.join(dir_sub, path_in_sub)
                        cases.add_case(Case(path_case))
                        break

print(f"{len(cases)} cases found:")
for case in cases:
    print("\n   " + case.name)
    print("      DC-bias:", case.dc_bias)
    print("      custom waveforms:", case.contains_custom)
    print("      phase:", case.cwaveform_phase)


xy = cases.get_value_pair("cwaveform_phase", "dc_bias", custom_waveform_only=True)
print(xy)

# generate figure folder if non existent
path_figures = os.path.join(dir_root, "Figures")
if not os.path.isdir(path_figures):
    os.mkdir(path_figures)

plotting.set_plot_globals()
plotting.plot_dc_bias_over_phase(xy)

# create XTplots
for case in cases:
    movie2xt(case.name, case.path_movie1, path_figures, case.freq*1e-6, num_zones=case.rffac+1)
