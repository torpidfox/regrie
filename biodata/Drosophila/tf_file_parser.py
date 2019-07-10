import pandas as pd
pd.set_option('display.max_columns', None)

tfs = "cad hkb kni gt tll bcd hb Kr".split()
# copy_numbers = [0] * len(tfs)

# for i, tf in enumerate(tfs):
# 	if tf == 'kni':
# 		copy_numbers[i] = 10
# 	if tf == 'gt':
# 		copy_numbers[i] = 1


with open('tf_100.csv') as f:
	df = pd.read_csv(f, sep=';')


print(df)
print(df[['name', 'REPRESSIONPROBABILITY', 'REPLENRIGHT', 'COPYNUMBER', 'REPRESSOR']])


is_repressor = {'bcd' : 'false', 'cad' : 'false', 'Kr' : 'false', 'hb' : 'true', 'gt' : 'true', 'kni' : 'true'}
copy_number = {'Kr' : 4, 'gt' : 6, 'hb' : 19, 'cad' : 7, 'tll' : 6, 'kni' : 5, 'hkb' : 1, 'bcd' : 5}
#copy_number = {tf : 100 for tf in tfs}
#copy_number = dict(zip(tfs, copy_numbers))



for index, row in df.iterrows(): 
	#df.set_value(index, 'COPYNUMBER', 1)
	# if row['name'] in copy_number.keys():
	# 	df.set_value(index, 'COPYNUMBER', copy_number[row['name']])
	# else:
	# 	df.set_value(index, 'COPYNUMBER', 0)
	if row['name'] in is_repressor.keys():
		#df.set_value(index, 'COPYNUMBER', 1)
		df.set_value(index, 'REPRESSOR', is_repressor[row['name']])

		if is_repressor[row['name']] == 'true':
			df.set_value(index, 'REPRESSIONPROBABILITY', 1)
			df.set_value(index, 'REPLENRIGHT', 125)
			df.set_value(index, 'REPLENLEFT', 125)
		else:
			df.set_value(index, 'REPRESSIONPROBABILITY', 0)
			df.set_value(index, 'REPLENRIGHT', 0)
			df.set_value(index, 'REPLENLEFT', 0)

#print(df[['name', 'REPRESSIONPROBABILITY', 'REPLENRIGHT', 'COPYNUMBER', 'REPRESSOR']])

#df.to_csv(path_or_buf='tf_1.csv', sep=';', index=False)
