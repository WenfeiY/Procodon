"""
Compute cracking complexity 
"""
### Available codecs
#----------------------------------------------------------------------------------------------------------------
from Procodon.codec.codec import Codecgenerator

codecgenerate = Codecgenerator('S.cerevisiae.codec_0.json', 'S.cerevisiae.codon_frequency.csv', 'S.cerevisiae.aa_frequency_dict.json')

codecgenerate.get_and_write_filtered_codecs('S.cerevisiae.codec_pass_dict.json')

for cut_off in [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]:
    filtered_codecs = codecgenerate.get_filtered_codecs(cut_off)
    print(cut_off)
    print(len(filtered_codecs))
#----------------------------------------------------------------------------------------------------------------


### Compute amount of possible combinations with interference genes
#----------------------------------------------------------------------------------------------------------------
'''
1 dubious gene: y = x + 1
2 dubious gene: y = (x+2)*(x+1)/2
3 dubious gene: y = (x+3)*(x+2)*(x+1)/6
4 dubious gene: y = (x+4)*(x+3)*(x+2)*(x+1)/24
5 dubious gene: y = (x+5)*(x+4)*(x+3)*(x+2)*(x+1)/120
'''
import numpy as np

def combinatorial_function(x, num):
    if num == 1:
        return x + 1
    if num == 2:
        return (x+2)*(x+1)/2
    if num == 3:
        return (x+3)*(x+2)*(x+1)/6
    if num == 4:
        return (x+4)*(x+3)*(x+2)*(x+1)/24
    if num == 5:
        return (x+5)*(x+4)*(x+3)*(x+2)*(x+1)/120
    
x = np.linspace(5, 200, 2000)
x.sort()

for num in [1, 2, 3, 4, 5]:
    y = combinatorial_function(x, num)
    with open(str(num) + '_dubious_gene.txt', 'w') as hd:
        for i in range(2000):
            hd.write(str(x[i]) + '\t' + str(y[i]) + '\n')
#----------------------------------------------------------------------------------------------------------------


### Compute amount of possible orders
#----------------------------------------------------------------------------------------------------------------
import math
import pandas as pd

#  Generate plot data
x_data = list(range(1, 21))        # X from 1 to 20
y_data = [math.factorial(x) for x in x_data]  # Calculate Y = X!

#  Save as csv
df = pd.DataFrame({'X': x_data, 'Y': y_data})
df.to_csv('factorial_data.csv', index=False)
#----------------------------------------------------------------------------------------------------------------
