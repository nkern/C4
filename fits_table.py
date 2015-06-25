import astropy.io.fits as fits
import numpy as np


def fits_table(dictionary,keys,filename):
	''' Takes data as a dictionary and makes a fits table, keys is a list or array 
		with the column names in the order you want them to be added to the fits file
		data in dictionary should be a numpy.ndarray '''
	length = len(keys)
	cols = []
	for i in range(length):
		format = dictionary[keys[i]].dtype.name
		if 'int' in format:
			form = 'J'
		elif 'float' in format:
			form = 'D'
		elif 'string' in format:
			form = 'J'
		elif 'bool' in format:
			form = 'L'
		cols.append(fits.Column(name=keys[i],format=form,array=dictionary[keys[i]]))
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto(filename)
	return

def fits_data(fits_data):
	d = {}
	for i in range(len(fits_data.dtype)):
		d[fits_data.columns[i].name] = fits_data[fits_data.columns[i].name][np.where(fits_data[fits_data.columns[i].name]!=0)]
	return d



