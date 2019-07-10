import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import re
from collections import namedtuple, defaultdict
from itertools import chain
from scipy.stats.stats import pearsonr
from textwrap import wrap
import operator
import mplcursors
from scipy.stats import ttest_ind, ks_2samp, wilcoxon
from sklearn.mixture import GaussianMixture
import pprint
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVC
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

pp = pprint.PrettyPrinter(indent=4)
#from statsmodels.stats import multitest

avail_bnd = [4582, 5108, 
8005, 8052, 
8113, 8778, 
8808, 8865, 
8900, 8937, 
8958, 9040,
9071, 9262, 
9309, 9492,
9508, 9734,
9775, 9850,
9871, 9905,
9930, 10218,
10275, 10398,
10420, 10502,
10564, 11999,
15724, 15752,
16786, 16805,
17916, 18000]

count = range(50, 1900, 50)
time = 120

tfs = "hb Kr cad bcd gt kni hkb tll".split()
is_repressor = {'bcd' : False, 'cad' : False, 'Kr' : False, 'hb' : True, 'gt' : True, 'kni' : True}
sites_count = {'bcd' : 31, 'cad' : 48, 'Kr' : 25, 'hb' : 135, 'gt' : 38, 'kni' : 31}
ignore_names = ['tll', 'hkb']
rep_len = 125

path_thres100 = '/run/media/alisa/Elements/test/results/server/results_thres_repression_1'
path_nothres100 = '/run/media/alisa/Elements/test/results/server/results_no_thres_repression_1/'

path_thres1 = '/run/media/alisa/Elements/test/results/server/results_thres_repression_1'
path_nothres1 = '/run/media/alisa/Elements/test/results/server/results_no_thres_repression_1/'

Ensemble = namedtuple('Ensemble', ['site', 'first_reached', 'times_reached', 'time_occupied', 'never_reached'])
Site = namedtuple('Site', ['tf', 'left', 'right', 'dir'])

def border_distance(site):
	i = 0

	while site.left > avail_bnd[i + 1]:
		i += 2

	#return avail_bnd[i + 1] - avail_bnd[i]
	return min(site.left - avail_bnd[i], avail_bnd[i + 1] - site.right)


def size_of_open_region(site):
	i = 0

	while site.left > avail_bnd[i + 1]:
		i += 2

	return avail_bnd[i + 1] - avail_bnd[i]

def is_repressed(site, sites):
	count = 0

	for neighbour in sites:
		if neighbour != site and is_repressor[neighbour.tf]:
			if site.left > neighbour.left - rep_len and site.right < neighbour.right + rep_len:
				count += 1

	return count


def overlap_count(site, sites):
	count = 0
	overlap = 0
	for neighbour in sites:
		if neighbour != site and (site.tf != neighbour.tf or\
		 (site.tf == neighbour.tf and site.dir != neighbour.dir)):
			if (site.right >= neighbour.right and site.left <= neighbour.right):
				count += 1
			elif (site.left <= neighbour.left and site.right >= neighbour.left):
				count += 1

	return count


def is_overlapping(site, sites):
	count = 0
	overlap = 0
	for neighbour in sites:
		if neighbour != site and (site.tf != neighbour.tf or\
		 (site.tf == neighbour.tf and site.dir != neighbour.dir)):
			if (site.right >= neighbour.right and site.left <= neighbour.right) and \
			(neighbour.right - site.left) > overlap:
				overlap = neighbour.right - site.left
			elif (site.left <= neighbour.left and site.right >= neighbour.left) and \
			(site.right - neighbour.left) > overlap:
				overlap = site.right - neighbour.left

	return overlap


def find_closest(site, sites):
	min_dist = float('inf')
	closest_site = None

	for neighbour in sites:
		if neighbour != site and neighbour.tf != site.tf:
			dist =  abs(neighbour.left - site.left)

				# if (site.right < neighbour.left and site.left > neighbour.left)\
				#  or (site.left > neighbour.left and site.left < neighbour.right):
				# 	return float('inf')

			if dist < min_dist and dist > (site.right - site.left):
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
	unvisited_count = 0

	with open(filename) as f:
		#skip header
		f.readline()

		for l in f:
			nums = l.split(',')

			if nums[4] == ' true\n' and nums[0] != '"hb:chr2R:5098..5108:0"':
				names.append(parse_site(nums[0].replace('"', '')))
				time_reached = float(nums[1])

				if time_reached == -1:
					unvisited_count += 1

				time_occupied = float(nums[3])
				count_reached = float(nums[2])

				time_reached, time_occupied = time_reached, time_occupied 
				reached_stats.append(time_reached)
				occupied_stats.append(time_occupied)
				times_reached.append(count_reached)
				never_reached.append(False)

	return [Ensemble(n, reached, reached_count, occ, was_reached) \
	for n, reached, reached_count, occ, was_reached in zip(names, reached_stats, times_reached, occupied_stats, never_reached)],\
	unvisited_count


def get_data(path, time=None):
	data = list()
	unvisited_counts = []
	print(time)
	for root, dirs, files in os.walk(path):
		if root.find('set') != -1:
			for f in files:
				if not time:
					if re.findall(r'target_site_\d', f) and '.0s' not in f:
						# with open(root + '/' + f) as file:
						# 	df = pd.read_csv(file)
						# data.append(df)
						d, unvisited_count = parse_file(root + '/' + f)
						data += d
						unvisited_counts.append(unvisited_count)
				else:
					if 'target_site_{}.0s'.format(time) in f:
						d, unvisited_count = parse_file(root + '/' + f)
						data += d
						unvisited_counts.append(unvisited_count)

	return data, unvisited_counts


def time_series_data(path):
	time = count

	data = []

	for t in time:
		d, _ = get_data(path, t)
		occ, t1, t2, _ = concat(d)
		occ, _, _ = get_mean(occ, t1, t2)
		data.append(occ)

	return data

def delta(t1, t2, site, time):
	return t1[site] / time - t2[site] / time


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
	unvisited_count = {x : 0 for x in names}

	for row in data:
		if row.first_reached != -1 and row.first_reached < 5000:
			means_reached[row.site].append(row.first_reached)
		if row.first_reached == -1:
			unvisited_count[row.site] += 1

		means_occupied[row.site].append(row.time_occupied)
		means_times_reached[row.site].append(row.times_reached)

	return means_occupied, means_reached, means_times_reached, unvisited_count

def get_mean(occupied, reached, times_reached):
	means_occupied = {name: np.mean(x) for name, x in occupied.items()}
	means_reached = {name: np.mean(x) for name, x  in reached.items()}
	means_times_reached = {name: np.mean(x) for name, x  in times_reached.items()}

	return means_occupied, means_reached, means_times_reached

def draw(data, label, color='r'):
	sns.kdeplot(data, color=color)
	plt.hist(data, label=label, density=True, alpha = 0.5, color=color)
	plt.xlabel('Время первого достижения, с',
		fontsize=16)
	plt.ylabel('Частота', 
		fontsize=16)

def sdo_ado(occupancy, energies, reached, condition='', names=None):
	if not names:
		names = [name for name in occupancy.keys() if condition in name.tf and name.right != 5108]

	print(energies)

	ado = [energies[name] for name in occupancy.keys() if name in names]

	reached_stats = [reached[name] for name in occupancy.keys() if name in names]

	sdo = [occupancy[name] for name in occupancy.keys() if name in names]
	
	ado /= max(np.abs(ado))

	ado = np.log(ado)
	#sdo = occupancy
	#print(list(zip(reached_stats, sdo, ado)))
	sdo /= max(sdo)
	#sdo /= max(sdo)
	sdo = np.log(sdo)
	sdo_vs_ado = [(x, y) for name, x, y in zip(names, sdo, ado) if name in names]

	# sdo_visited_a_lot = [(x, y) for times_reached, x, y in zip(reached_stats, sdo, ado) if times_reached > thres]
	# sdo_visited_less = [(x, y) for times_reached, x, y in zip(reached_stats, sdo, ado) if times_reached < thres]
	# names_visited_less = [name for name, times_reached in reached.items() if times_reached > thres and condition in name.tf]
	# print('Overvisited')
	# print(names_visited_less)

	return sdo_vs_ado, names

def draw_scatterplot(ado, sdo, xlabel='Вероятность занятости по формуле', 
	ylabel='Вероятность занятости из симуляций', 
	label=None, title=None, col='r', names=None):

	plt.title('Соотношение времен первого достижения')
	to_draw = [(x, y) for x, y in zip(ado, sdo) if x + 100 > y]
	ado, sdo = zip(*to_draw)
	plt.scatter(ado, sdo, label=label, alpha=0.5)
	plt.xlabel(xlabel, fontsize=16)
	plt.ylabel(ylabel, fontsize=16)

	# mplcursors.cursor().connect("add", 
	# 	lambda sel: sel.annotation.set_text(names[sel.target.index]))

	# if names:
	# 	for i, txt in enumerate(names):
	# 		plt.annotate(txt, (ado[i], sdo[i]))



	borders = [min(min(sdo), min(ado)),
	 max(max(ado), max(sdo))]
	plt.plot(borders, borders, color='k')
		

def stat_test(d1, d2, test_mixture=False):
	print(ttest_ind(d1, d2, equal_var=False))
	#print(ks_2samp(d1, d2))
	print(wilcoxon(d1, d2))

	if test_mixture:
		mixture1 = GaussianMixture(n_components=2).fit([[el] for el in d1])
		print(mixture1.means_)
		mixture2 = GaussianMixture(n_components=2).fit([[el] for el in d2])
		print(mixture2.means_)

def count_outliers(stat):
	count = {tf : 0 for tf in tfs}
	out = list()

	for site, time in stat.items():
		if time > 1300:
			count[site.tf] += 1
			out.append(site)

	for tf in tfs:
		time_mean = [val for el, val in stat.items() if el.tf == tf]
		print(min(time_mean))
		count[tf] /= len(time_mean)
		print(tf, np.mean(time_mean), len(time_mean))

	print(count)
	pp.pprint(out)

def repression_test():
	d, unvisited_count_no_repression = get_data(path_nothres)
	d, unvisited_count_repression = get_data(path_thres)

	results, edges = np.histogram(unvisited_count_no_repression, density=True)
	binWidth = edges[1] - edges[0]
	plt.bar(edges[:-1], results*binWidth, binWidth,
		alpha=0.5,
		color='r',
		label='Без репрессии')

	#sns.kdeplot(thres[s], color='g')
	results, edges = np.histogram(unvisited_count_repression, density=True)
	binWidth = edges[1] - edges[0]
	plt.bar(edges[:-1], results*binWidth, binWidth,
		alpha=0.5,
		color='g',
		label='C репрессией')

	plt.legend()
	plt.title('Распределение количества непосещенных сайтов')

	plt.show()

	print(ttest_ind(unvisited_count_no_repression,
	 unvisited_count_repression,
	  equal_var=False))

	print(np.mean(unvisited_count_no_repression),
		np.mean(unvisited_count_repression))



def classify(x, y):
	x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2)
	reg = SVC().fit(x_train, y_train)
	print(reg.score(x_test, y_test))

def plot_relative_differ():
	d, _ = get_data(path_nothres100)
	occ1, reached1, times_reached, _ = concat(d)
	occupied_stats, reached_stats_no_thres100, times_reached_stats = get_mean(occ1, reached1, times_reached)

	d, _ = get_data(path_thres100)
	occ1, reached1, times_reached, _ = concat(d)
	occupied_stats, reached_stats_thres100, times_reached_stats = get_mean(occ1, reached1, times_reached)

	d, _ = get_data(path_nothres1)
	occ1, reached1, times_reached, _ = concat(d)
	occupied_stats, reached_stats_no_thres1, times_reached_stats = get_mean(occ1, reached1, times_reached)

	d, _ = get_data(path_thres1)
	occ1, reached1, times_reached, _ = concat(d)
	occupied_stats, reached_stats_thres1, times_reached_stats = get_mean(occ1, reached1, times_reached)
	names = [name for name in reached_stats_thres100.keys() if name.tf not in ignore_names]

	reached_stats_no_thres100 = [reached_stats_no_thres100[name] for name in names]
	reached_stats_thres100 = [reached_stats_thres100[name] for name in names]
	reached_stats_no_thres1 = [reached_stats_no_thres1[name] for name in names]
	reached_stats_thres1 = [reached_stats_thres1[name] for name in names]

	deltas100, deltas1 = [], []
	for q1, q2 in zip(reached_stats_thres100, reached_stats_no_thres100):
		deltas100.append((q2 - q1) / q2)

	for q1, q2 in zip(reached_stats_thres1, reached_stats_no_thres1):
		deltas1.append((q2 - q1) / q2) 

	plt.scatter(deltas1, deltas100)
	plt.xlabel('Относительная разница в случае 8 молекул')
	plt.ylabel('Относительная разница в случае 600 молекул')
	plt.show()

def first_reached_hist():
	d, _ = get_data(path_nothres100)
	occ1, reached1, times_reached, _ = concat(d)
	occupied_stats, reached_stats, times_reached_stats = get_mean(occ1, reached1, times_reached)
	names = [name for name in reached_stats.keys()]

	#log_reached1 = np.log([reached_stats[name] for name in names])
	log_reached1 = [reached_stats[name] for name in names]


	draw(log_reached1, 
		label='Исходная модель')
	#print(np.array2string(log_reached1, separator=','))
	#print(reached_stats)

	d, _ = get_data(path_thres100)
	occ2, reached2, times_reached2, _ = concat(d)

	count = {tf : 0 for tf in tfs}



	# for site in reached1.keys():
	# 	r = ttest_ind(reached1[site], reached2[site], equal_var=False)

	# 	if r.statistic > r.pvalue:
	# 		print(site, 'passed')
	# 	else:
	# 		count.update({site.tf:count[site.tf] + 1})
	# 		print(site, 'did not pass')

	#print(count)


	occupied_stats, reached_stats2, times_reached_stats = get_mean(occ2, reached2, times_reached2)
	count_outliers(reached_stats)
	log_reached2 = [reached_stats2[name] for name in names]

	draw(log_reached2, 
		label='Модифицированная модель', color='g')

	plt.legend(fontsize=16)
	plt.title('\n'.join(wrap('Распределение логарифмов времен первого достижения, 100 молекул каждого белка')))
	plt.show()

	tfs_to_draw = [tf for tf in tfs]

	#draw_scatterplot(log_reached1, log_reached2, xlabel='Без порога', ylabel='С порогом')
	reg = LinearRegression().fit([[el] for el in log_reached1], log_reached2)
	xs = np.arange(min(min(log_reached1), min(log_reached2)), max(max(log_reached1), max(log_reached2)), 0.1)
	ys = [reg.coef_ * x + reg.intercept_ for x in xs]
	# plt.plot(xs, ys, label='Регрессионная прямая', color='k',
	# 	linestyle='--')


	for tf in tfs_to_draw:
		to_draw = [(el1, el2) for name, el1, el2 in zip(names, log_reached1, log_reached2) if name.tf == tf]
		xs, ys = list(zip(*to_draw))

		draw_scatterplot(xs, ys, xlabel='Время в исходной модели, с', 
			ylabel='Время в модифицированной модели, с',
			label=tf)

	borders = [min(min(log_reached1), min(log_reached2)),
	 max(max(log_reached1), max(log_reached2))]
	plt.plot(borders, borders, color='k')

	# energies = parse_energy('/run/media/alisa/Elements/test/reGRiE/pwm_energies.txt')
	# energies_filtered = [score for name, score in energies.items() if name in names]

	overlap_flag = list(map(lambda s: is_overlapping(s, names), names))
	overlap_counts = list(map(lambda s: overlap_count(s, names), names))
	#repress_count = list(map(lambda s: is_repressed(s, names), names))

	print('weirdos')
	for name, r1, r2, c in zip(names, log_reached1, log_reached2, overlap_flag):
		if r1 < r2:
			print(name, c)

	#x = list(zip(overlap_flag, repress_count, overlap_counts))
	#classify(x, [1 if r1 > r2 else 0 for r1, r2 in zip(log_reached1, log_reached2)])

	# non_overlapping_sites1 = [t for t, flag in zip(log_reached1, repress_count) if flag < 1]
	# non_overlapping_sites2 = [t for t, flag in zip(log_reached2, repress_count) if flag < 1]
	# print(len(non_overlapping_sites1))

	# draw_scatterplot(non_overlapping_sites1, non_overlapping_sites2, 
	# 	xlabel='Исходная модель', 
	# 	ylabel='Модифицированная модель',
	# 	label='Сайты без пересечений')

	# non_overlapping_sites1 = [t for t, flag in zip(log_reached1, repress_count) if flag >= 1]
	# non_overlapping_sites2 = [t for t, flag in zip(log_reached2, repress_count) if flag >= 1]

	# draw_scatterplot(non_overlapping_sites1, non_overlapping_sites2, 
	# 	xlabel='Исходная модель', 
	# 	ylabel='Модифицированная модель',
	# 	label='Сайты с пересечениями')

	#print(np.array2string(log_reached2, separator=','))

	#print([name for name, val in reached_stats2.items() if val > np.exp(7)])

	print('means distribution')
	stat_test(log_reached1, log_reached2, test_mixture=True)
	plt.legend(fontsize=12)
	#print(np.mean(np.exp(log_reached1)))
	#print(np.mean(np.exp(log_reached2)))

	plt.show()
	
	#s1 = sorted(zip(log_reached1, names), key=lambda x: x[1].left)
	#s2 = sorted(zip(log_reached2, names), key=lambda x:x[1].left)

	s1 = sorted(zip(log_reached1, names), key=lambda x: x[1].left)
	s2 = sorted(zip(log_reached2, names), key=lambda x: x[1].left)
	s1 = log_reached1
	s2 = log_reached2

	#xs = [x.left for _, x in s1]

	deltas = []
	to_ignore = []

	for d1, d2 in zip(s1, s2):
		q1 = d1
		q2 = d2

		deltas.append((q1 - q2))

	s = sorted(zip(deltas, s1), key=lambda x: x[1])
	d_s, fr_s = zip(*s)
	means = np.convolve(d_s, np.ones((10,))/10, mode='valid')

	# sl1 = [el for el in s1 if el > b1]
	# sl2 = [el for el in s2 if el > b2]

	# quantiles1 = [np.quantile(sl1, x) for x in [0.5, 1]]
	# quantiles2 = [np.quantile(sl2, x) for x in [0.5, 1]]

	# for q1, q2 in zip(quantiles1, quantiles2):
	# 	a1 = [el for el in sl1 if el < q1]
	# 	a2 = [el for el in sl2 if el < q2]		

	# 	deltas.append(np.mean(np.asarray(a1) - np.mean(a2)) / np.mean(a1))


	# sl1 = [el for el in s1 if el > borders1[-1]]
	# sl2 = [el for el in s2 if el > borders2[-1]]

	# deltas.append(np.mean(np.asarray(sl1) - np.mean(sl2)) / np.mean(sl1))

	plt.scatter(fr_s, d_s)
	#plt.scatter(sl1[:293], deltas[:293], label='Доля разницы между моделями')
	#plt.plot(fr_s[:len(means)], means, label='Скользящее среднее для N=10', color='r')
	plt.legend()
	plt.xlabel('Время первого достижения в исходной модели', fontsize=16)
	plt.ylabel('Разница между моделями', fontsize=16)
	plt.title('1 молекула')


	
	#plt.savefig('../logs/first_reached_100_repression.png', figsize=(20, 100), dpi=160)

def delta_plot():
	time_series_nonthres = time_series_data(path_nothres)
	time_series_thres = time_series_data(path_thres)
	deltas = []
	print(time_series_thres)
	print(time_series_nonthres)

	for i, t in enumerate(count):
		print(t)
		deltas.append([delta(time_series_thres[i], time_series_nonthres[i], s, t) for s in time_series_thres[i].keys()\
		 if s in time_series_nonthres[i].keys()])


	delta_disp = [np.mean(d) for d in deltas]
	
	#plt.hist(deltas[-1])
	#plt.show()
	#plt.clf()
	print(delta_disp)

	plt.plot(list(count), delta_disp)
	plt.xlabel('Время, с')
	plt.ylabel(r'$\delta$')
	plt.title('Влияние порога на времена занятости')
	plt.show()

def occ_hist(thres, nonthres, site, t):
	s = parse_site(site)
	#sns.kdeplot(nonthres[s], color='r')
	means = [np.mean(x) for name, x in nonthres.items() if name.tf not in ignore_names]
	results, edges = np.histogram(means, density=True)
	binWidth = edges[1] - edges[0]
	plt.bar(edges[:-1], results*binWidth, binWidth,
		alpha=0.5,
		color='r',
		label='Исходная модель')

	#sns.kdeplot(thres[s], color='g')
	means = [np.mean(x) for name, x in thres.items() if name.tf not in ignore_names]
	results, edges = np.histogram(means, density=True)
	binWidth = edges[1] - edges[0]
	plt.bar(edges[:-1], results*binWidth, binWidth,
		alpha=0.5,
		color='g',
		label='Модифицированная модель')
	
	plt.title('\n'.join(wrap('Распределение времен занятости сайта {} для времени {} с'.format(site, t))),
		fontsize=20)
	plt.legend(fontsize=14)
	plt.xlabel('Время, с', fontsize=16)
	plt.ylabel('Частота', fontsize=16)
	#plt.savefig('../logs/occ_{}_at_{}.png'.format(s.tf, t), figsize=(100, 500), dpi=160)
	plt.show()
	plt.clf()

def count_neighbour(site, sites, r=40):
	count = 0
	prev_site = site
	distance = 0

	for neighbour in sites:
		if neighbour != site and neighbour.tf == site.tf:
			if site.right > neighbour.left and site.left > neighbour.left or\
			site.right < neighbour.left and site.right < neighbour.left:
				#if neighbour.right - prev_site.right > 1:
				if not ((prev_site.right >= neighbour.right and prev_site.left <= neighbour.right) or\
					(prev_site.left <= neighbour.left and prev_site.right >= neighbour.left)):
					if abs(site.left - neighbour.left) <= r or abs(neighbour.right - site.right) <=r:
						count += 1
						prev_site = neighbour

	return count


	#return min(site.left - avail_bnd[i], avail_bnd[i + 1] - site.right)

def find_closest_in_open_region(site, sites):
	op1 = size_of_open_region(site)
	count = 0

	for s in sites:
		if s.tf == site.tf and size_of_open_region(s) == op1 and s != site:
			count += 1

	return count

#def find_farest



def unvisited_count_vs_overlap():
	d, _ = get_data(path_nothres)
	_, _, _, visited_count = concat(d)

	neighbour_count = list(map(lambda x: size_of_open_region(x), visited_count.keys()))
	#neighbour_count = list(map(size_of_open_region, visited_count.keys()))

	xs, ys = zip(*sorted(zip(neighbour_count, visited_count.values()), key=lambda x: x[0]))
	mean_x = {el : list() for el in set(xs)}

	for x, y in zip(xs, ys):
		mean_x[x].append(y)

	for k in mean_x.keys():
		mean_x.update({k : np.mean(mean_x[k])})

	xs, ys = zip(*sorted(zip(mean_x.keys(), mean_x.values()), key=lambda x: x[0]))
	print(xs)
	plt.step(xs, [y / 1000 for y in ys]) 
	plt.xlabel('Размер (в п. о.) области открытого хроматина')
	plt.ylabel('Доля симуляций, в которых сайт не был посещен')
	plt.title('\n'.join(wrap('Зависимость частоты достигаемости сайтов от доступности области открытого хроматина')))
	plt.show()


def first_reached_vs_bnd_distance():
	d, _ = get_data(path_thres)
	occ1, reached1, times_reached, _ = concat(d)
	occupied_stats, reached_stats, times_reached_stats = get_mean(occ1, reached1, times_reached)

	sites = [s for s in reached_stats.keys()]
	vals = [val for name, val in reached_stats.items()]

	neighbour_count = list(map(lambda x: is_overlapping(x, sites), sites))

	xs, ys = zip(*sorted(zip(neighbour_count, vals), key=lambda x: x[0]))
	print(ys)
	mean_x = {el : list() for el in set(xs)}

	for x, y in zip(xs, ys):
		mean_x[x].append(y)

	for k in mean_x.keys():
		mean_x.update({k : np.mean(mean_x[k])})

	xs, ys = zip(*sorted(zip(mean_x.keys(), mean_x.values()), key=lambda x: x[0]))
	print(xs)
	print(ys)
	plt.plot(xs, [y for y in ys]) 
	plt.xlabel('Количество пересечений с сайтами других белков')
	plt.ylabel('Время первого достижения',
		fontsize=16)
	plt.title('\n'.join(wrap('Зависимость времени достижения от взимного расположения сайтов')))
	plt.show()



def first_reached_vs_overlap():
	d, _ = get_data(path_thres1)
	occ1, reached1, times_reached, _ = concat(d)
	occupied_stats, reached_stats, times_reached_stats = get_mean(occ1, reached1, times_reached)

	#log_reached1 = np.log([el for name, el in reached_stats.items()])

	neighbour_count = defaultdict(None)
	sites_count = defaultdict()
	tfs_mean = [tf for tf in tfs]

	for tf in tfs_mean:
		sites = list(filter(lambda x: tf in x.tf, reached_stats.keys()))
		neighbours = list(map(lambda x: count_neighbour(x, sites), sites))
		print(tf, neighbours)
		n = [n for n in neighbours if n > 100]

		neighbour_count[tf] = np.mean(neighbours)
		#neighbour_count[tf] = np.mean(neighbours)

	#b = [1914.9989696296027, 922, 1777, 1883, 1000, 835, 718, 1265]
	#b = [922, 1777, 1883, 1000, 835, 1265]
	print(neighbour_count)

	xs, ys, names = zip(*sorted(zip(neighbour_count.values(), b, tfs_mean), key=lambda x: x[0]))

	for i, txt in enumerate(names):
		plt.annotate(txt, (xs[i], ys[i]), fontsize=18)

	plt.plot(xs, ys)
	#plt.title('Зависимость времени первого достижения от скученности сайтов')
	plt.xlabel('Среднее количество соседей в радиусе 40 п. о.',
		fontsize=16)
	plt.ylabel('Среднее время первого достижения',
		fontsize=16)
	plt.show()
	# no_hb = [el for el in reached_stats.keys() if el.tf != 'kni'']
	# neighbour_count_no_hb = list(map(lambda x: count_neighbour(x, no_hb), no_hb))
	# print(np.mean(neighbour_count_no_hb))

	# no_hb = [el for el in reached_stats.keys() if el.tf == 'kni']
	# neighbour_count_no_hb = list(map(lambda x: count_neighbour(x, no_hb), no_hb))
	# print(np.mean(neighbour_count_no_hb))

	# xs, ys = zip(*sorted(zip(neighbour_count, log_reached1), key=lambda x: x[0]))
	# mean_x = {el : list() for el in set(xs)}

	# for x, y in zip(xs, ys):
	# 	mean_x[x].append(y)

	# for k in mean_x.keys():
	# 	mean_x.update({k : np.mean(mean_x[k])})

	# plt.plot(mean_x.keys(), mean_x.values()) 
	# plt.show()


a = [135, 25, 48, 31, 38, 30, 7, 37]
b = [1914.9989696296027, 922, 1777, 1883, 1000, 835, 718, 1265]
xs, ys = zip(*sorted(zip(a, b), key=lambda x: x[0]))
#plt.plot(xs, ys)
#plt.show()



#plot_relative_differ()
#delta_plot()
#first_reached_vs_overlap()
#first_reached_vs_bnd_distance()
first_reached_hist()
#unvisited_count_vs_overlap()
#a = [3.5, 1.5, 1.9 , 1.3, 0.7, 0.8,0.07,0.8]
#plt.plot(xs, ys)
#plt.show()

# d, _ = get_data(path_thres)
# occ1, reached1, times_reached, _ = concat(d)
# occupied_stats, reached_stats, times_reached_stats = get_mean(occ1, reached1, times_reached)
# energies = parse_energy('/run/media/alisa/Elements/test/reGRiE/pwm_energies.txt')
# sites = [el for el in energies.keys() if el.tf not in ignore_names and el.right != 5108]

#repression_test()
#unvisited_count_vs_overlap()
#first_reached_hist()
#delta_plot()
# occ1, _, _, _ = concat(get_data(path_nothres100, 600)[0])
# occ2, _, _, _ = concat(get_data(path_thres100, 600)[0])
# occ_hist(occ2, occ1, 'Kr:chr2R:9623..9634:1', 1000)

#first_reached_hist()


# for i, tf in enumerate(tfs):
# 	if tf not in ignore_names:
# 		sdo_vs_ado, names = sdo_ado(energies, occupied_stats, times_reached_stats, condition=tf)
# 		# overlap_flag = list(map(lambda s: overlap_count(s, energies.keys()), names))
# 		# to_draw_1 = [el for el, overlap in zip(sdo_vs_ado, overlap_flag) if overlap == 0]
# 		# if i == 0:
# 		# 	draw_scatterplot(*zip(*to_draw_1), col='b', label='Сайты без пересечений')
# 		# else:
# 		# 	draw_scatterplot(*zip(*to_draw_1), col='b')

# 		to_draw_1 = [el for el in sdo_vs_ado if el[0] < el[1] + 0.3]
# 		# if i == 0:
# 		# 	draw_scatterplot(*zip(*to_draw_1), col='g', label='Сайты c пересечениями')
# 		# else:
# 		draw_scatterplot(*zip(*to_draw_1), col='g')

# overlap_flag = list(map(lambda s: overlap_count(s, energies.keys()), energies.keys()))
# non_overlapping_sites = [site for site, flag in zip(energies.keys(), overlap_flag) if flag == 0 and site in sites]
# sdo_vs_ado, names = sdo_ado(occupied_stats, energies,times_reached_stats, names=non_overlapping_sites)
# draw_scatterplot(*zip(*sdo_vs_ado), label='Сайты, не имеющие пересечений с другми', names=non_overlapping_sites, col='b')

# overlapping_sites = [site for site, flag in zip(energies.keys(), overlap_flag) if flag > 0 and site in sites]
# sdo_vs_ado, names = sdo_ado(occupied_stats, energies, times_reached_stats)
# draw_scatterplot(*zip(*sdo_vs_ado), label='Сайты, имеющие пересесения с другими', col='r')


# quick = [name for name, time in reached_stats.items() if np.log(time) < 7]
# print(quick)

# to_draw = [el for el, name in zip(sdo_vs_ado, energies.keys()) if name in overlapping_sites and site in sites]
# draw_scatterplot(*zip(*to_draw), label='сайты, имеющие пересечения с другими сайтами')
# plt.legend()
# plt.show()

#overlap_hist(parse_site('gt:chr2R:9984..9996:1'), parse_site('kni:chr2R:9990..10003:1'), reached)


#sdo, _ = zip(*sdo_a_lot)
#print(list(zip(energies.keys(), sdo, distances)))
#print(weirdos_nothres)
#plt.hist(list(times_reached_stats.values()))
#draw(np.log([el for el in reached_stats.values() if el != -1]), label='Без репрессии')
#draw_occupancy(energies, occ)
#print(weirdos_nothres)


# for site in reached.keys():
# 	if len([el for el in reached[site] if el != -1]) > 1 and len([el for el in reached2[site] if el != -1]) > 1:
# 		draw(np.log([el for el in reached[site] if el != -1]), label='Без порога', color='r')
# 		draw(np.log([el for el in reached2[site] if el != -1]), label='С порогом', color='g')
# 		plt.legend()
# 		plt.show()


# for tf in tfs:
# 	sdo_vs_ado= sdo_ado(energies, occupied_stats, times_reached_stats, condition=tf)
# 	draw_scatterplot(*zip(*sdo_vs_ado),  title='', label=tf)
#sdo2, ado, weirods_thres = sdo_ado(energies, occupied_stats)
#draw_scatterplot(sdo2, ado)

# d1 = [0.33197120446422423, 0.40100342661903726, 0.4039449004141973, 0.4138218569728289, 0.4064170707321282, 0.37060797549293223, 0.37945326340070307, 0.36877717147756534, 0.36982836535072916, 0.37038494558984947, 0.3771872645281038, 0.37081432045685686, 0.3701028447002683, 0.3730916061479226, 0.38360347368085207, 0.37896870529627075, 0.3713099235355674, 0.37089904242772637, 0.37340016710622564, 0.41691365386851503, 0.4233041653450488, 0.4062590018185312, 0.4261856631986548, 0.4181704992603954, 0.4281284915145199, 0.43034558896399056, 0.4236416647470908, 0.422517162059087, 0.4256179864679118, 0.4217460424950903, 0.415824634503199, 0.41479910317892854, 0.41001739966533013, 0.41227149996111007, 0.41085024723354974, 0.4072005857888806, 0.40717981507721573, 0.3905431141298288, 0.39285731803591256, 0.3838189842502899, 0.38924349986259477, 0.3768398187680393, 0.38265520432657113, 0.37001890678380034, 0.3810865910528344, 0.396487763935799, 0.3960235644346286, 0.39691624188858265, 0.39677313591164315, 0.39467516786999407, 0.3951003562205149, 0.39926627649113605, 0.3882594796932211, 0.3951651684592587, 0.39446082993963477, 0.4019619428615835, 0.40982925131267434, 0.4042313247060931, 0.40386754720033713, 0.41424157613504725, 0.41837405222954854, 0.4187605217256465, 0.417519771703553, 0.42645223064596016, 0.4196833278226427, 0.4229498073263219, 0.4144990705822802, 0.43589753233331424, 0.43215871242975906, 0.44146629239022783, 0.44015430749076667, 0.43883295077950474, 0.40572143776784175, 0.3969920307615179, 0.41720622710002725, 0.41580495760373387, 0.4148147089020964, 0.4220129983870441, 0.42778182265387216, 0.39455876849627836, 0.3956641779335315, 0.387791939417449, 0.3925217973254788, 0.4001529846846615, 0.39274367447903374, 0.39397887446604146, 0.38971697636873076, 0.3972601609847372, 0.3928216554891221, 0.3766927698342863, 0.3656815945163699, 0.35410278534352035, 0.355917281042385, 0.35047926357067505, 0.34701145654536875, 0.32354004825764876, 0.3192407830020423, 0.3147763702970119, 0.29660640445391473, 0.30216176228678243, 0.2967331907084137, 0.29585410888784786, 0.30684394026709155, 0.3171065056458112, 0.3070969391205013, 0.30301106718672827, 0.30973100235868306, 0.2932922178312992, 0.29836142177592573, 0.28876097257728545, 0.27351425018410347, 0.27658090650069284, 0.271815573105185, 0.26373452183780594, 0.2628434948198874, 0.26455129289261675, 0.26122973488354917, 0.2613036419961436, 0.2516959466910903, 0.23327474899553877, 0.232858625474819, 0.2317465514184351, 0.24136045515784182, 0.22907694272624515, 0.23395213903688752, 0.2303211281816004, 0.22116213972552895, 0.21469816457855892, 0.2072296782063598, 0.20198247387107898, 0.20033013132932553, 0.20291561739109007, 0.20265360145133104, 0.20936724795212533, 0.26038626957035926, 0.24911641162584153, 0.24957623308380544, 0.24869897931518734, 0.259562353450422, 0.24366225423027807, 0.24432123354168953, 0.2450286639051916, 0.22270788748684328, 0.22182171173357554, 0.2229007113410783, 0.21956850497945157, 0.22264114025740442, 0.21603261205878127, 0.20004845932998594, 0.20063239930070786, 0.2009033757189386, 0.2151315154857302, 0.21328681647538594, 0.23001578485273916, 0.2316180796852229, 0.23524675026181266, 0.24267443736821534, 0.2296693039820834, 0.2006739705680639, 0.21173078959092906, 0.21148274718053833, 0.2157650130958797, 0.2314875829843425, 0.18278519932399065, 0.21225203064217862, 0.20066603947923725, 0.20312371196515167, 0.2161026252129427, 0.21992151011399175, 0.241845355303793, 0.2279291187525525, 0.23034093758088786, 0.22614605566186854, 0.22896574631373234, 0.2241561149668896, 0.20615083924694377, 0.19832530074684096, 0.2054449795342199, 0.20284051566940703, 0.19957685002166847, 0.2093627866986513, 0.21377896814676067, 0.20549207597076477, 0.19371186070495572, 0.1879759694778141, 0.18784855559136754, 0.1638205299479055, 0.16569665421713156, 0.15776652292825985, 0.15234228463647004, 0.1272333647942803, 0.13476066056997224, 0.1166795188647463, 0.10005232602198458, 0.11060673156579061, 0.10372573651690782, 0.07560761557841636, 0.07209472112303679, 0.09258530935206362, 0.09062687558983691, 0.08816828141427338, 0.0861466622237768, 0.08296628434046595, 0.06115372917667751, 0.05788038134856832, 0.08123912080913732, 0.08757265805084578, 0.08927710757409384, 0.07823338664693213, 0.07602042626783596, 0.0575342462404445, 0.03650402161247795, 0.039073621998374275, 0.0383593776802973, 0.04019264762954202, 0.030929312556240894, 0.031221382868618187, 0.0397960471477044, 0.0458006560349699, 0.04550302484284208, 0.04063553595903119, 0.05945410772503921, 0.05706194896122615, 0.06728591535430394, 0.08121209968403997, 0.07050135686443404, 0.06328032496201037, 0.06566499196264604, 0.01287616053273878, 0.004792508434043215, 0.014051359935368416, 0.00010137675372217174, 0.027129691355934064, 0.0285247594722688, 0.029438122708708696, 0.027289828669633524, 0.027379507159590935, 0.018595814936469404, 0.01459976396917685, 0.012355675537311128, 0.01055682540121514, 0.00982712736550055, 0.0005287521598181704, 0.015764762876195087, 0.010298136179446212, 0.01311535808572379, 0.011671604169570662, 0.022265054869541064, 0.018177146035013074, 0.0223069661565457, 0.030015903848416157, 0.023110051702298386, 0.02465457658440445, 0.02599795268549179, 0.02492960895636734, 0.03388404020382063, 0.04079754640718164, 0.04279513569587438, 0.04403992906921208, 0.04738887354666635, 0.052476869166337745, 0.047057156648302186, 0.048285877444148716, 0.04873192118081924, 0.05121000650071398, 0.04980405004124727, 0.050382911685646604, 0.052432909265104674, 0.055927569684300016, 0.05459972477708571, 0.05966123686846078, 0.06691362216501194, 0.06720199542497012, 0.06305577310540919, 0.06119363739567487, 0.06568248312239656, 0.0677436471901231, 0.06094607563575099, 0.0641500926950483, 0.05620784145542863, 0.05451006449077498, 0.037724078643212734, 0.0219554415082638, 0.03652246100185915, 0.03230733781221906, 0.05345220260582139, 0.057088769923923255, 0.05557081541616856, 0.0659722800928115, 0.06225876720455153, 0.06798419891188144, 0.07285238836312391, 0.06757068810636821, 0.06915910666353815, 0.06906713872603705, 0.08258275461189711, 0.08173871784989213, 0.07633858696795297, 0.06410867811347677, 0.03770992332562757, 0.038670680477081915, 0.012431557451612995, 0.007411460949741863, 0.023713588668887148, 0.17033047281905, 0.17086195103814086]
# d2 = [0.040651151856613894, 0.016083936505355, 0.02874662823007658, 0.0339731421085598, 0.035766685815243865, 0.035564369669139645, 0.03971792903951541, 0.04030706450983141, 0.03825874031973457, 0.033364844531516294, 0.03143930285492654, 0.022588322781929, 0.021433790848620467, 0.028670972449984844, 0.041837449386704334, 0.043206413920229184, 0.029709000699737864, 0.017738020608969552, 0.021262610517859937, 0.018885083562801722, 0.020921147858613844, 0.023043367780109437, 0.023247536074363723, 0.021481513589438724, 0.02412286698170977, 0.0204630428423784, 0.02378299537157078, 0.019861390361867043, 0.019964615777500375, 0.019992497497380628, 0.01954566237530973, 0.020578781232913422, 0.019905499461394485, 0.023244522462239983, 0.021993840006678686, 0.023168376301486177, 0.019745473250934995, 0.028889556171661987, 0.0262341411422818, 0.023597576259132582, 0.01902849582859248, 0.01753757323733509, 0.01430661926430391, 0.013961965995178667, 0.016869344651861233, 0.018866260932228575, 0.017152717407761767, 0.017241195481205947, 0.036474895669403505, 0.03921779526627035, 0.037163228076874305, 0.04169263638956556, 0.04475998557342775, 0.040219922756049964, 0.032987514028486764, 0.030998232578986547, 0.033492474505514486, 0.03396041784625471, 0.035180066316171, 0.0410591231424055, 0.04323632042392701, 0.014875778594058375, 0.0003215870279245783, 0.0007497064961415249, 6.706945011226269e-05, 0.002846989271790496, 0.009976400680606899, 0.020135225887130938, 0.025688273818728517, 0.02962828171022083, 0.03364868887151622, 0.021679207923941297, 0.02344428336131722, 0.025611712652733938, 0.08172629234354412, 0.12660507179788907, 0.02162545997466808, 0.02597343438419514, 0.009831600779255342, 0.005597098084447941, 0.0005471582028708983, 0.06866165306811851, 0.08142514833106124, 0.07714580444549782, 0.028350726111719942, 0.022679909997538418, 0.021990876801317318, 0.011262673485701308, 0.08514095578780707, 0.07867227292338803, 0.07643653975594536, 0.06637311636974018, 0.05868010710971026, 0.049343667606799924, 0.04877385282928132, 0.04625379302591997, 0.046750597630814204, 0.04155654203557191, 0.041032875687419534, 0.04328088322148885, 0.0433668852697565, 0.042526031843080674, 0.04150864402388373, 0.040599477376824525, 0.043030618950273235, 0.044105290787698696, 0.042236915933019176, 0.04130395526657977, 0.04326023879299839, 0.04019386791214873, 0.03953423196576863, 0.03676647637008558, 0.038009360325214306, 0.039834096977235126, 0.037685599456824447, 0.038775336068001647, 0.040616616691496216, 0.044650466393012754, 0.04710816299231681, 0.046961758949134104, 0.04745068569756029, 0.05012610089423153, 0.05573079285297214, 0.05586306290300782, 0.05249709724843942, 0.051565921306670315, 0.05149598833195367, 0.05248356072906404, 0.04920632555191677, 0.04664506573664101, 0.045054771844628604, 0.04505691430948997, 0.04528375248665116, 0.04273553756375485, 0.044728142719394855, 0.04412242157489015, 0.043818598899208354, 0.04217848729815977, 0.04111791521864484, 0.040909647752467006, 0.04116367267664058, 0.040242598573975, 0.04016788927115518, 0.03964349780928538, 0.03944508528616331, 0.039639946799887665, 0.04010982127797186, 0.03853051446782764, 0.03781348754604229, 0.037361231848538935, 0.03613268866005494, 0.036375140282660645, 0.037532637035844335, 0.03622862716030329, 0.03634768680545467, 0.03616355607393106, 0.03587338764826998, 0.03493507448953861, 0.03506804667992002, 0.035008174238428794, 0.03161230960687097, 0.028862849419632887, 0.02600866918921758, 0.02535789673550705, 0.025258365107325462, 0.02509850307757439, 0.025369418055461293, 0.0242151576336348, 0.02352235620222648, 0.02308880947242862, 0.024121597121303805, 0.023781907091134438, 0.02315814211137277, 0.023389519026925917, 0.024407168296409383, 0.024525685142157524, 0.023841960024353896, 0.023914899272669685, 0.023157002593824805, 0.023278081701281903, 0.023156484220478734, 0.023273844663581213, 0.022768904329790125, 0.021996106636707403, 0.022325528488652067, 0.02287074015403321, 0.022443367450236444, 0.021818975990909793, 0.020897419993677054, 0.020160705101267187, 0.019864518022107674, 0.02078267655498225, 0.020449606355171825, 0.019245524192339592, 0.019148584015291643, 0.01950873294802733, 0.020474516110508826, 0.020346523266799164, 0.018705614083677328, 0.0186693635363833, 0.019155473652620915, 0.019241980957212225, 0.01893166611326796, 0.019076029134072615, 0.01853714985686366, 0.018207340153800607, 0.01843200562599866, 0.01681095848052303, 0.016076380701298843, 0.015867536674306917, 0.016007272309849427, 0.016056273958940676, 0.01562562614475095, 0.014510724151123415, 0.014431059968846508, 0.01415138730778356, 0.013870171350595506, 0.015397052316354179, 0.01567864143032905, 0.01591460695131912, 0.015789827827516686, 0.01592746008921976, 0.0162052440336711, 0.016889299226722172, 0.017380815285184353, 0.015680178565049892, 0.015360830274397847, 0.015052404897661515, 0.01472721095737206, 0.014741908921036174, 0.015283039106594374, 0.015135428289074177, 0.014779194457693741, 0.015244853768361762, 0.013869147241245935, 0.014538303369731102, 0.014645686458346087, 0.014367990361407823, 0.01596752728090827, 0.015046133892269812, 0.015189462553416518, 0.013427995871895656, 0.014133212660708043, 0.012213728445051452, 0.011922837111538782, 0.012209978325281649, 0.010475707181352746, 0.012261172166791807, 0.011639876810014504, 0.010226183391360565, 0.010884053928989328, 0.012034879289407231, 0.012942357658040233, 0.014417055430820055, 0.015140095349642994, 0.01604171647276584, 0.017191143514317116, 0.019230363275673047, 0.02064713401078752, 0.020577295877161905, 0.022340318424659954, 0.01871031938490853, 0.017621441147023853, 0.016439912812445265, 0.017835296759975226, 0.01750689192913301, 0.015630549755121747, 0.016585779259532075, 0.019810721192735418, 0.02070299172592437, 0.020549490338314592, 0.020936544213144094, 0.02161169464062775, 0.023631983964895742, 0.025543002549918335, 0.024844543712873493, 0.018853227884182564, 0.01969097112769553, 0.01886197954406528, 0.02573646330042357, 0.020391261077895222, 0.025352772626568466, 0.03025312129565761, 0.028923773388490535, 0.02819159615946891, 0.030039602168533224, 0.029399620635034106, 0.026595877930268064, 0.029280791243414164, 0.03360048487055112, 0.03355749248639611, 0.0412020625666929, 0.04146595661982386, 0.049838402771111304, 0.0472392422812032, 0.04466201414197962, 0.04296173302906415, 0.03922457152609867, 0.04648756553651123, 0.04502903887670868, 0.03057125709409102, 0.0343964617560077] 
# print(len(d1))
# print(len(d2))
# plt.scatter(d2, d1[:302])
# plt.xlabel('Относительная разница в случае 8 молекул')
# plt.ylabel('Относительная разница в случае 600 молекул')
# plt.show()

plt.legend()
plt.show()


