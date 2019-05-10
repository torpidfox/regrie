import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import re
from collections import namedtuple, defaultdict
from itertools import chain
from scipy.stats.stats import pearsonr
from textwrap import wrap


avail_bnd = [4582, 5108, 8005, 8052, 8113, 8778, 8808, 8865, 8900, 8937, 8958, 9040,
 9071, 9262, 9309, 9492, 9508, 9734, 9775, 9850, 9871, 9905, 9930, 10218, 10275, 10398, 10420,
  10502, 10564, 11999, 15724, 15752, 16786, 16805, 17916, 18000]

time = 120

tfs = "cad hkb kni gt tll bcd hb".split()

path_nothres = '/run/media/alisa/Elements/test/reGRiE/results_debug_no_thres_no_repression_100/'
path_thres = '/run/media/alisa/Elements/test/reGRiE/results_debug_no_thres_repression_100/'
Ensemble = namedtuple('Ensemble', ['site', 'first_reached', 'times_reached', 'time_occupied', 'never_reached'])
Site = namedtuple('Site', ['tf', 'left', 'right', 'dir'])

def is_overlapping(site, sites):
	for neighbour in sites:
		if neighbour != site:
			if (site.right > neighbour.left and site.left < neighbour.left)\
				 or (site.left > neighbour.left and site.left < neighbour.right):
					return True

	return False

def find_closest(site, sites):
	min_dist = float('inf')
	closest_site = None

	for neighbour in sites:
		if neighbour != site:
			if neighbour.dir == site.dir and neighbour.tf != site.tf:
				dist =  abs(neighbour.left - site.left)

				# if (site.right < neighbour.left and site.left > neighbour.left)\
				#  or (site.left > neighbour.left and site.left < neighbour.right):
				# 	return float('inf')

				if dist < min_dist and dist < (site.right - site.left):
					min_dist = dist
					closest_site = neighbour


	return min_dist#, closest_site


def parse_site(site):
	chunks = site.split(':')
	tf = chunks[0]
	left, right = chunks[2].split("..")
	direction = chunks[-1]

	return Site(tf, int(left), int(right), int(direction))


def parse_file(filename, log=True):
	names, reached_stats, times_reached, occupied_stats, never_reached = list(), list(), list(), list(), list()

	with open(filename) as f:
		#skip header
		f.readline()

		for l in f:
			nums = l.split(',')

			if nums[4] == ' true\n':
				names.append(parse_site(nums[0].replace('"', '')))
				time_reached = float(nums[1])
				time_occupied = float(nums[3])
				count_reached = float(nums[2])

				time_reached, time_occupied = time_reached, time_occupied 
				reached_stats.append(time_reached)
				occupied_stats.append(time_occupied)
				times_reached.append(count_reached)
				never_reached.append(False)

	print(len(list(filter(lambda x: x != None, never_reached))))

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


def plot_heatmap(sdo, reached_stats, distance):
	reached_stats_values = np.asarray(list(reached_stats.values()))
	print(sdo)
	#print(np.concatenate(([np.asarray(sdo)], [reached_stats_values], [np.asarray(distance)]), axis=0))
	sns.heatmap(np.concatenate(([sdo], [reached_stats_values], [np.asarray(distance)]), axis=1))


def parse_energy(filename):
	points, names = list(), list()

	with open(filename) as f:
		for l in f:
			row = l.split()
			names.append(parse_site(row[0]))
			points.append(float(row[-1]))

	ys = list(map(lambda x: 15.0 * np.exp(x), points))

	return {name : e for name, e in zip(names, ys)}

def draw_energy(ys):
	ax = sns.kdeplot(np.log(ys), color='b')
	ax.set_label('Expected occupancy times distribution')
	plt.hist(np.log(ys), label='Expected times distribution', density=True, alpha = 0.5, color='b')

def concat(data):
	names = [x.site for x in data]
	means_reached = {x : list() for x in names}
	means_occupied = {x : list() for x in names}
	means_times_reached = {x : list() for x in names}

	for row in data:
		means_reached[row.site].append(row.first_reached)
		means_occupied[row.site].append(row.time_occupied)
		means_times_reached[row.site].append(row.times_reached)

	return means_occupied, means_reached, means_times_reached

def get_mean(occupied, reached, times_reached):
	means_occupied = {name: np.mean(x) for name, x in occupied.items()}
	means_reached = {name: np.mean(x) for name, x  in reached.items()}
	means_times_reached = {name: np.mean(x) for name, x  in times_reached.items()}

	return means_occupied, means_reached, means_times_reached

def draw(data, label, color='r'):
	sns.kdeplot(data, color=color)
	plt.hist(data, label=label, density=True, alpha = 0.5, color=color)

def sdo_ado(energies, occupancy, reached, condition=''):
	names = [name for name in occupancy.keys() if condition in name]
	ado = [energies[name] for name in occupancy.keys() if condition in name.tf]
	reached_stats = [reached[name] for name in occupancy.keys() if condition in name.tf]

	sdo = [occupancy[name] for name in occupancy.keys() if condition in name.tf]
	

	ado /= max(np.abs(ado))

	ado = np.log(ado)
	#sdo = occupancy
	#print(list(zip(reached_stats, sdo, ado)))
	sdo /= max(sdo)
	#sdo /= max(sdo)
	sdo = np.log(sdo)

	weirdos = [(name, s, a) for name, s, a in zip(occupancy.keys(), sdo, ado) if a > s]
	print(weirdos)

	sdo_vs_ado = [(x, y) for name, x, y in zip(names, sdo, ado) if condition in name]

	# sdo_visited_a_lot = [(x, y) for times_reached, x, y in zip(reached_stats, sdo, ado) if times_reached > thres]
	# sdo_visited_less = [(x, y) for times_reached, x, y in zip(reached_stats, sdo, ado) if times_reached < thres]
	# names_visited_less = [name for name, times_reached in reached.items() if times_reached > thres and condition in name.tf]
	# print('Overvisited')
	# print(names_visited_less)

	return sdo_vs_ado
def draw_scatterplot(sdo, ado, xlabel='ADO', ylabel='SDO', label=None, title=None, col='r'):
	plt.title('ADO vs SDO, 100 molecules of each tf' + title)
	plt.scatter(ado, sdo, label=label)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.plot( [-7,1], [-7,1])

def draw_occupancy(energies, occupancy, weirdos):
	bnd = [parse_bnd(name) for name in occupancy.keys()]

	ys = list(chain.from_iterable([(occ, occ) for occ in occupancy.values()]))
	e = list(chain.from_iterable([(occ, occ) for occ in energies.values()]))
	xs = list(chain.from_iterable([(left, right) for left, right in bnd]))
	site_names = list(chain.from_iterable([(name, name) for name in energies.keys()]))

	sort = sorted(zip(xs, ys, e, site_names), key=lambda x: x[0])
	xs, ys, e, names = list(zip(*sort))
	ys = np.asarray([x[0] for x in ys])
	ys /= max(np.abs(ys))
	e /= max(np.abs(e))

	weirdos_ys = [y for name, y in zip(names, ys) if name in weirdos]
	weirdos_xs = [x for name, x in zip(names, xs) if name in weirdos]

	plt.scatter(xs, ys, label='Результат симуляции')
	plt.scatter(xs, e, label='ADO')
	plt.scatter(weirdos_xs, weirdos_ys)

def reached_from_distance(distances, reached):
	reached_vs_distance = defaultdict(list)

	for d, s in zip(distances, reached):
		reached_vs_distance[d].append(s)

	for k, v in reached_vs_distance.items():
		reached_vs_distance.update({k : np.mean(v)})

	to_plot = sorted(reached_vs_distance.items(), key=lambda x: x[0])
	x, y = zip(*to_plot)
	print(to_plot)

	plt.plot(x, y)

def overlap_hist(site1, site2, times_reached):
	frac = [np.log(s1 / s2) for s1, s2 in zip(times_reached[site1], times_reached[site2]) if s1 != -1 and s2 != -1]
	sns.kdeplot(frac)



plotType = 'Total occupancy'
plotType = 'First reached'

names = ['Распределение времен занятости', 'Распределение времен первого достижения']
x_names = ['Логарифм времени занятости', 'Натуральный логарифм времени первого достижения']

energies = parse_energy('/run/media/alisa/Elements/test/reGRiE/pwm_energies.txt')

#draw_energy(list(energies.values()))
#plt.title(names[0] if plotType == 'Total occupancy' else names[1])
#plt.xlabel(x_names[0] if plotType == 'Total occupancy' else x_names[1])
#plt.ylabel('Частота')
#plt.legend()
#plt.savefig('/run/media/alisa/Elements/dying/d/logs/total_occupancy_mean_log.png')


occ1, reached, times_reached = concat(get_data(path_nothres))
occupied_stats, reached_stats, times_reached_stats = get_mean(occ1, reached, times_reached)
#draw(np.log([el if el != -1 else 1e7 for el in reached_stats.values()]), label='Без порога')

# for tf in tfs:
# 	sdo_vs_ado= sdo_ado(energies, occupied_stats, times_reached_stats, condition=tf)
# 	draw_scatterplot(*zip(*sdo_vs_ado), title=', no thres', col='g')
# overlap_flag = list(map(lambda s: is_overlapping(s, energies.keys()), energies.keys()))
# overlapping_sites = [site for site, flag in zip(energies.keys(), overlap_flag) if flag]
#print(overlapping_sites)

#overlap_hist(parse_site('gt:chr2R:9984..9996:1'), parse_site('kni:chr2R:9990..10003:1'), reached)


#sdo, _ = zip(*sdo_a_lot)
#print(list(zip(energies.keys(), sdo, distances)))
#print(sorted(distances))
#reached_from_distance(*list(zip(*occupied_crossing)))
#plot_heatmap(sdo, times_reached_stats, distances)
#print(weirdos_nothres)
#plt.hist(list(times_reached_stats.values()))
draw(np.log([el for el in reached_stats.values() if el != -1]), label='Без репрессии')
#draw_occupancy(energies, occ)
#print(weirdos_nothres)

#draw(np.log([el for el in reached[parse_site('gt:chr2R:5022..5034:1')] if el != -1]), label='Без порога', color='r')
occ2, reached, times_reached = concat(get_data(path_thres))
#draw(np.log([el for el in reached[parse_site('gt:chr2R:5022..5034:1')] if el != -1]), label='С порогом', color='g')
occupied_stats, reached_stats, times_reached_stats = get_mean(occ2, reached, times_reached)


# for tf in tfs:
# 	sdo_vs_ado= sdo_ado(energies, occupied_stats, times_reached_stats, condition=tf)
# 	draw_scatterplot(*zip(*sdo_vs_ado),  title='', label=tf)
#sdo2, ado, weirods_thres = sdo_ado(energies, occupied_stats)
#draw_scatterplot(sdo2, ado)
#print(weirods_thres)

draw(np.log([el if el != -1 else 1e7 for el in reached_stats.values()]), label='С репрессией', color='g')
#draw_occupancy(energies, occ1, weirods_thres)
#sdo, ado = sdo_ado(energies, occ1)

#draw_scatterplot(sdo1, sdo2, xlabel='С порогом', ylabel='Без порога')

#plot_heatmap(sdo)


#reached_from_distance(distances, times_reached_stats.values())
#reached_vs_bnd_dist(times_reached_stats)
plt.legend()
plt.title('Распределение времен первого достижения')
#plt.show()

plt.savefig('../logs/first_reached_repression_100mol.png', figsize=(20, 100), dpi=160)

for i in range(0, len(avail_bnd), 2):
	print(avail_bnd[i], avail_bnd[i+1])