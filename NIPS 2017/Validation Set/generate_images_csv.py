import glob, os
import numpy as np
import csv

os.chdir('images/')

try:
    # initialize
    folders = [name for name in os.listdir(".") if os.path.isdir(name)]
    image_ids      = np.zeros(0, dtype=str)
    true_labels    = np.zeros(0, dtype=int)
    true_label = 0
    
    # record image ids and true labels while moving images out
    for folder in folders:
        true_label += 1
        os.chdir(folder + '/')
        try:
            for file in glob.glob("*.JPEG"):
                image_ids =   np.append(image_ids,   file[:-5])
                true_labels = np.append(true_labels, true_label)
                os.rename(file, "../" + file)
        finally:
            os.chdir('../')
        # os.rmdir(folder + '/')

finally:
    os.chdir('../')

# randomly assign target class and remove collision
target_classes = np.random.randint(1, 1001, size=np.size(true_labels))
for i in np.where(target_classes == true_labels)[0]:
    while target_classes[i] == true_labels[i]:
        target_classes[i] = np.random.randint(1, 1001, size=1)

# save image info
with open('images.csv', 'w') as csvfile:
    writer = csv.DictWriter(csvfile, ['ImageId','URL','x1','y1','x2','y2','TrueLabel',
    								  'TargetClass','OriginalLandingURL','License',
    								  'Author','AuthorProfileURL'],
    						lineterminator='\n')
    writer.writeheader()
    for image_id, true_label, target_class in \
    zip(image_ids, true_labels, target_classes):
        writer.writerow({'ImageId': image_id,
                         'TrueLabel': true_label,
                         'TargetClass': target_class})