"""This script should be called by "diffdef.sh". It operates on
different lines with common variables. It requires version 2.7 or
higher of Python (it also works with Python 3).

Author: Lionel GUEZ"""

import sys, os.path

with open(os.path.join(sys.argv[1], "common_var")) as common_var1, \
        open(os.path.join(sys.argv[2], "common_var")) as common_var2, \
        open("common_lines", "a") as common_lines, \
        open("uniq_val", "w") as uniq_val:
    uniq_val.write("Comparison of " + sys.argv[1] + " and " +  sys.argv[2] 
                   + ":\n") # title line
    for line1, line2 in zip(common_var1, common_var2):
        i = line1.index("=")
        var_name = line1[:i]
        value1 = line1[i + 1:-1]
        value2 = line2[i + 1:-1]
        try:
            value1_num = float(value1)
            value2_num = float(value2)
        except ValueError:
            uniq_val.write(var_name + "=\t" + value1 + "\t" + value2 + "\n")
        else:
            if (value1_num, value2_num) == (0, 0) \
                    or abs(value2_num - value1_num) \
                    / max(abs(value1_num), abs(value2_num)) <= 1e-5:
                # These should not be integer values, write a formatted value:
                common_lines.write(var_name + "=" + str(value1_num) + "\n")
            else:
                # These could be integer values, write the original strings:
                uniq_val.write(var_name + "=\t" + value1 + "\t" + value2 
                               + "\n")
