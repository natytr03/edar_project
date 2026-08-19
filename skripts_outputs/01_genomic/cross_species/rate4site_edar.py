import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

df = pd.read_csv(
    r"/data/results/r4s_output.res",
    comment='#',
    sep=r'\s+',
    header=None,
    names=['pos', 'aa', 'score', 'qq_low', 'qq_high', 'std', 'msa']
)

#smoothing line shows general direction of the data
df['smooth'] = df['score'].rolling(window=5, center=True).mean()

fig = plt.figure(figsize=(20, 6))
gs = gridspec.GridSpec(2, 1, height_ratios=[5, 0.6], hspace=0.05)

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharex=ax1)

ax1.plot(df['pos'], df['score'], alpha=0.3, linewidth=0.8, color='steelblue')
ax1.plot(df['pos'], df['smooth'], linewidth=2, color='steelblue')
ax1.axvline(x=370, color='red', linewidth=1.5, label='V370A')
ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)

ax1.set_ylabel('Evolutionary rate (Rate4Site)')
ax1.set_title('EDAR conservation profile (Rate4Site)')
ax1.legend(loc='upper left', fontsize=8)
ax1.tick_params(labelbottom=False)

#show the conserved block human as reference, domain boundaries from UniProt/InterPro
domains = [
    ("CRD1",         30, 72,   "#a6cee3"),
    ("CRD2",         73, 114,  "#1f78b4"),
    ("CRD3",         115, 149, "#b2df8a"),
    ("TM",           185, 209, "#fdbf6f"),
    ("Death Domain", 344, 436, "#fb9a99"),
]

for label, start, end, color in domains:
    ax2.barh(y=0, width=end-start, left=start, height=0.8, color=color, edgecolor='black', linewidth=0.5)
    ax2.text((start+end)/2, 0, label, ha='center', va='center', fontsize=7, rotation=0)

ax2.axvline(x=370, color='red', linewidth=1.5)  #highlight V370A
ax2.set_ylim(-0.5, 0.5)
ax2.set_yticks([])
ax2.set_xlabel('Protein position')

plt.savefig('rate4site_profile_domains_track.png', dpi=300, bbox_inches='tight')
plt.show()