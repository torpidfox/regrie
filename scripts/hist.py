import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import re
from collections import namedtuple

avail_bnd = [4582, 5108, 8005, 8052, 8113, 8778, 8808, 8865, 8900, 8937, 8958, 9040,
 9071, 9262, 9309, 9492, 9508, 9734, 9775, 9850, 9871, 9905, 9930, 10218, 10275, 10398, 10420,
  10502, 10564, 11999, 15724, 15752, 16786, 16805, 17916, 18000]

time = 120

path_nothres = '/run/media/alisa/Elements/test/reGRiE/results_debug/'
path_thres = '/run/media/alisa/Elements/test/reGRiE/results_debug_thres_30/'
Ensemble = namedtuple('Ensemble', ['name', 'first_reached', 'times_reached', 'time_occupied', 'never_reached'])


def parse_file(filename, log=True):
	names, reached_stats, times_reached, occupied_stats, never_reached = list(), list(), list(), list(), list()

	with open(filename) as f:
		#skip header
		f.readline()

		for l in f:
			nums = l.split(',')

			#if nums[0][1:3] == 'Kr':

			if nums[4] == ' true\n':
				names.append(nums[0].replace('"', ''))
				time_reached = float(nums[1])
				time_occupied = float(nums[3])
				count_reached = float(nums[2])

				if log and time_occupied != 0:
					time_reached, time_occupied = np.log(time_reached), time_occupied 
					reached_stats.append(time_reached)
					occupied_stats.append(time_occupied)
					times_reached.append(count_reached)
					never_reached.append(False)
				elif time_occupied == 0:
					#occupied_stats.append(float('NaN'))
					never_reached.append(True)

	print(len(list(filter(lambda x: x, never_reached))))

	return [Ensemble(n, reached, reached_count, occ, was_reached) \
	for n, reached, reached_count, occ, was_reached in zip(names, reached_stats, times_reached, occupied_stats, never_reached)]


def get_data(path):
	data = list()
	for root, dirs, files in os.walk(path):
		if root.find('set') != -1:
			for f in files:
				if re.findall(r'target_site_\d+', f):
					# with open(root + '/' + f) as file:
					# 	df = pd.read_csv(file)
					# data.append(df)
					data += parse_file(root + '/' + f)

	return data

def parse_bnd(description):
	bnd = description.split(':')[2]
	left, _, right = bnd.split('.')

	return int(left), int(right)

def is_available(description):
	left, right = parse_bnd(description)

	for i in range(0, len(avail_bnd), 2):
		if right < avail_bnd[i + 1] and left > avail_bnd[i]:
			return True

	return False

def parse_energy(filename):
	points, names = list(), list()

	with open(filename) as f:
		for l in f:
			row = l.split()
			names.append(row[0])
			#if row[0] == 'Kr':
			points.append(float(row[-1]))


	ys = list(map(lambda x: 15 * np.exp(x), points))

	return {name : e for name, e in zip(names, ys)}

def draw_energy(ys):
	ax = sns.kdeplot(np.log(ys), color='b')
	ax.set_label('Expected occupancy times distribution')
	plt.hist(np.log(ys), label='Expected times distribution', density=True, alpha = 0.5, color='b')

def concat(data):
	names = set((x.name for x in data))
	means_reached = {x : list() for x in names}
	means_occupied = {x : list() for x in names}
	means_times_reached = {x : list() for x in names}

	for row in data:
		means_reached[row.name].append(row.first_reached)
		means_occupied[row.name].append(row.time_occupied)
		means_times_reached[row.name].append(row.times_reached)

	return means_occupied, means_reached, means_times_reached

def get_mean(occupied, reached, times_reached):
	means_occupied = {name: np.mean(x) for name, x in occupied.items()}
	means_reached = {name: np.mean(x) for name, x  in reached.items()}
	means_times_reached = {name: np.mean(x) for name, x  in times_reached.items()}

	return means_occupied, means_reached, means_times_reached

def draw(data, label, plotType='total occupancy', color='r', log=False):
	sns.kdeplot(data, color=color)
	plt.hist(data, label=label, density=True, alpha = 0.5, color=color)
	#ax.set_title(plotType)
	#plt.xlim([0, 150])

def draw_scatterplot(energies, occupancy):
	ado = [energies[name] for name in occupancy.keys()]
	ado /= max(np.abs(ado))
	ado = np.log(ado)
	sdo = list(occupancy.values())
	sdo /= max(np.abs(sdo))
	sdo = np.log(sdo)


	plt.title('ADO vs SDO')
	plt.scatter(ado, sdo)
	plt.xlabel('ADO')
	plt.ylabel('SDO')
	plt.plot( [min(min(ado), min(sdo)),1], [min(min(ado), min(sdo)),1])



#plotType = 'Total occupancy'
plotType = 'First reached'

names = ['Распределение времен занятости', 'Распределение времен первого достижения']
x_names = ['Логарифм времени занятости', 'Натуральный логарифм времени первого достижения']

energies = parse_energy('/run/media/alisa/Elements/test/reGRiE/pwm_energies.txt')
draw_energy(list(energies.values()))
plt.title(names[0] if plotType == 'Total occupancy' else names[1])
plt.xlabel(x_names[0] if plotType == 'Total occupancy' else x_names[1])
plt.ylabel('Частота')
plt.legend()
#plt.savefig('/run/media/alisa/Elements/dying/d/logs/total_occupancy_mean_log.png')


occ, reached, times_reached = concat(get_data(path_nothres))
occupied_stats, reached_stats, times_reached_stats = get_mean(occ, reached, times_reached)
#draw_scatterplot(energies, occupied_stats)
#plt.hist(list(times_reached_stats.values()))
#draw(np.log(list(reached_stats.values())), label='Без порога')

occ, reached, times_reached = concat(get_data(path_thres))
occupied_stats, reached_stats, _ = get_mean(occ, reached, times_reached)
#draw_scatterplot(energies, occupied_stats)
draw(np.log(list(occupied_stats.values())), label='С порогом', color='g')
plt.legend()
plt.show()
