Calculation of the RMSF:

```python
import numpy as np
import MDAnalysis as mda

u = mda.Universe("topol.tpr", "trj.xtc")
ca = u.select_atoms("name CA")
means = np.zeros((len(ca), 3))
sumsq = np.zeros_like(means)
for k, ts in enumerate(u.trajectory):
    sumsq += (k/(k+1.0)) * (ca.positions - means)**2
    means[:] = (k*means + ca.positions)/(k+1.0)
rmsf = np.sqrt(sumsq.sum(axis=1)/(k+1.0))

matplotlib.pyplot.plot(ca.residues.resids, rmsf)
```

And plot more nicely....
```python
import matplotlib.pyplot as plt
import seaborn.apionly as sns
%matplotlib inline

plt.style.use('ggplot')
sns.set_style('ticks')

fig = plt.figure(figsize=(4,2))
ax = fig.add_subplot(111)
color = sns.color_palette()[2]
ax.fill_between(ca.residues.resids, rmsf, alpha=0.3, color=color)
ax.plot(ca.residues.resids, rmsf, lw=2, color=color)
sns.despine(ax=ax)
ax.set_xlabel("Residue")
ax.set_ylabel(r"C$_\alpha$ RMSF ($\AA$)")
ax.set_xlim(1, max(ca.residues.resids))
fig.savefig("ca_rmsf.pdf")
```