import os

class Case:
    def __init__(self):
        self.powerICP = 0
        self.powerCCP1 = 0
        self.powerCCP2 = 0


cases = []
newCase = Case()
cases.append(newCase)

print(cases[0].powerICP)

# root directory
dir_root = os.path.abspath("D:\\UIGEL5_D_Florian\Voltage_Waveform_Tailoring\HPEM\ArCF4O2\DarkSpaceGeometry\ConstVoltage_200_1000")
# dir_root = "D:\\UIGEL5_D_Florian\Voltage_Waveform_Tailoring\HPEM\ArCF4O2\DarkSpaceGeometry\ConstVoltage_200_1000"
dirs_in_root = os.listdir(dir_root)
print(dirs_in_root)

# print header

print("type", 10*" ", "name")
print(50*"-")
for dir_in_root in dirs_in_root:

    dir_sub = os.path.join(dir_root, dir_in_root)
    # check if file:
    string_dir_type = "{0:<17}".format("not recognised")
    if os.path.isfile(dir_sub):
        string_dir_type = "{0:<17}".format("file")

    # check if directory
    if os.path.isdir(dir_sub):
        string_dir_type = "{0:<17}".format("directory")

        # if directory, find out if source or case directory
        # ---------------------------------------------------
        paths_in_sub = os.listdir(dir_sub)
        for path_in_sub in paths_in_sub:
            #print(path_in_sub)
            if os.path.isfile(os.path.join(dir_sub, path_in_sub)):
                #print("here")
                #print(path_in_sub)

                if "icp.f" in path_in_sub.lower() or "icp.for" in path_in_sub.lower():
                    string_dir_type = "{0:<17}".format("source directory")

    # print(dir)
    print(string_dir_type, dir_in_root)

