import matplotlib.pyplot as plt
import squarify 

# 1. Define the Data
# Structure: {Category: {Disease: Value}}
data = {
    'Cardiovascular': {'Coronary heart disease': 40, 'Stroke': 15, 'Atrial': 10},
    'Mental disorders': {'Depressive disorders': 25, 'Anxiety disorders': 15, 'Bipolar': 8},
    'Neoplasms': {'Lung cancer': 20, 'Breast cancer': 12, 'Bowel cancer': 10, 'Pancreatic': 5},
    'Musculoskeletal': {'Back pain': 25, 'Osteoarthritis': 12},
    'Neurological': {'Dementia': 20, 'Migraine': 10, 'Epilepsy': 8},
    'Respiratory': {'COPD': 15},
    'Digestive': {'Chronic liver disease': 12},
    'Metabolic': {'Type 2 diabetes': 18},
    'Other': {'Hearing loss': 10}
}

# 2. Flatten data for squarify
labels = []
values = []
colors = []

# Color palette matching your image
category_colors = {
    'Cardiovascular': '#636EFA', 'Mental disorders': '#00CC96', 
    'Neoplasms': '#19D3F3', 'Musculoskeletal': '#FFA15A', 
    'Neurological': '#FF6692', 'Respiratory': '#FF97FF', 
    'Digestive': '#EF553B', 'Metabolic': '#AB63FA', 'Other': '#B6E880'
}

for category, sub_dict in data.items():
    for disease, val in sub_dict.items():
        labels.append(f"{disease}")
        values.append(val)
        colors.append(category_colors[category])

# 3. Plotting
plt.figure(figsize=(14, 6))
squarify.plot(sizes=values, label=labels, color=colors, alpha=0.8, 
              edgecolor="white", linewidth=2, text_kwargs={'fontsize':9, 'color':'#333'})

plt.title("Disease Burden Treemap", fontsize=16)
plt.axis('off')
plt.show()