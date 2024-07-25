import numpy
import config

well_labels = []
for letter in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
    for number in range(1,13):
        well_labels.append('%s%d' % (letter , number))


antibiotic_labels = ['cN',  'c11',  'c9',  'c7']
antibiotic_conc = 25*numpy.asarray([0, 2**-11, 2**-9,  2**-7])
print(antibiotic_conc)

antibiotic_labels_dict = {}
for a_idx, a in enumerate(antibiotic_labels):
    antibiotic_labels_dict[a] = antibiotic_conc[a_idx]

phage_labels = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7']
#phage_labels = phage_labels[::-1]
phage_dilution = numpy.asarray([10**(-1*int(s[1])) for s in phage_labels])
phage_labels_dict = {}
for a_idx, a in enumerate(phage_labels):
    phage_labels_dict[a] = phage_dilution[a_idx]



def parse_metadata():

    label_well_dict = {}

    line_count = 0
    for line in open('%slabels.csv' % config.data_directory,  encoding='utf-8-sig'):

        line = line.strip().split(',')

        for item_idx, item in enumerate(line):
            #label_well_dict[well_labels[(12*line_count) + item_idx]] = item
            label_well_dict[well_labels[(12*line_count) + item_idx]] = item

        line_count += 1

    return label_well_dict




def parse_data():

    data_dict = {}
    data_dict['hours'] = []

    label_well_dict = parse_metadata()

    for key, value in label_well_dict.items():
        data_dict[value] = []

    plate_data = open('%skitp_qbio_turner_2024_07_24_abx_phage_conc_gradient.csv' % config.data_directory, 'r',  encoding='utf-8-sig')
    header = plate_data.readline()
    header = header.strip().split(',')
    for line in plate_data:
        
        line = line.strip().split(',')
        line = [float(s) for s in line]
        data_dict['hours'].append(line[0]/(60*60))

        for od_idx, od in enumerate(line[1:]):
            
            data_dict[label_well_dict[header[od_idx+1]]].append(od)

    plate_data.close()

    for key, value in data_dict.items():

        data_dict[key] = numpy.asarray(value)

    return data_dict



#data_dict = parse_data()


#print(data_dict)

