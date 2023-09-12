import re

result = open('./output.txt', 'r')

test = result.read()

total_var_re = re.compile('Phased Vars: (.*?) variants')
overall_re = re.compile('Overall: (.*?) hap')
correct_re = re.compile('Overall.*total,(.*?)hap')

total_var = sum(list(map(lambda x: int(x.strip(' ').strip('\t')), total_var_re.findall(test))))
total_hap = sum(list(map(lambda x: int(x.strip(' ').strip('\t')), overall_re.findall(test))))
corre_hap = sum(list(map(lambda x: int(x.strip(' ').strip('\t')), correct_re.findall(test))))

print('Phased {} vars in total, resulting in {} haplotypes, {} of which is correct, with ACC: {}%'\
      .format(total_var, total_hap, corre_hap, round(corre_hap/total_hap*100, 4)))