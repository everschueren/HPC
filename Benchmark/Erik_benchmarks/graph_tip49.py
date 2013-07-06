import sys


data = open('AllScores.out','r')
D = data.readlines()
data.close()

saint,nsaf,zscore,dscore,mist = {},{},{},{},{}
for d in D:
    d = d.strip().split('\t')
    p = (d[0],d[1])
    saint[p] = float(d[5])
    nsaf[p] = float(d[6])
    zscore[p] = float(d[7])
    dscore[p] = float(d[8])
    mist[p] = float(d[9])

DB = [('UXT', 'LRPPRC'), ('NUFIP1', 'FMR1'), ('KIAA0515', 'EHMT2'), ('VPS72', 'EWSR1'), ('RP11-529I10.4', 'RUVBL2'), ('RUVBL2', 'RP11-529I10.4'), ('KIAA0515', 'EHMT2'), ('RUVBL2', 'PON2'), ('RUVBL1', 'RUVBL2'), ('RUVBL2', 'RUVBL1'), ('RUVBL1', 'C2orf44'), ('RUVBL2', 'C2orf44'), ('RUVBL1', 'RUVBL2'), ('RUVBL2', 'RUVBL1'), ('RUVBL2', 'USP7'), ('RUVBL2', 'RPS11'), ('RUVBL1', 'E2F7'), ('RP11-529I10.4', 'RUVBL1'), ('RUVBL1', 'RP11-529I10.4'), ('RP11-529I10.4', 'RUVBL2'), ('RUVBL2', 'RP11-529I10.4'), ('TP53', 'TP53BP1'), ('TP53', 'TP53BP1'), ('TP53', 'TP53BP1'), ('TP53', 'BRE'), ('TP53', 'BRE'), ('TP53', 'BRE'), ('TP53', 'BRCC3'), ('TP53', 'BRCC3'), ('TP53', 'TP53BP1'), ('TP53', 'TP53'), ('TP53', 'RPA1'), ('TP53', 'TP53BP1'), ('TP53', 'TP53BP1'), ('TP53', 'PRKDC'), ('TP53', 'SFN'), ('TP53', 'TP53BP1'), ('TP53', 'TP53'), ('TP53', 'RPA1'), ('TP53', 'BRCC3'), ('TP53', 'PRMT1'), ('TP53', 'YBX1'), ('TP53', 'PARP1'), ('TP53', 'XRCC1'), ('POLR2E', 'C19orf2'), ('TP53', 'RPL11'), ('TP53', 'RPL11'), ('TP53', 'MDC1'), ('TP53', 'MTA2'), ('TP53', 'EEF2'), ('NUFIP1', 'FMR1'), ('TP53', 'PRKDC'), ('RUVBL1', 'CTNNB1'), ('RUVBL1', 'CTNNB1'), ('RUVBL1', 'CTNNB1'), ('RUVBL1', 'CTNNB1'), ('RUVBL2', 'EP400'), ('RUVBL1', 'EP400'), ('RUVBL1', 'EP400'), ('RUVBL1', 'EP400'), ('RUVBL1', 'EP400'), ('TP53', 'HDAC1'), ('TP53', 'SMARCB1'), ('TP53', 'SMARCA4'), ('TP53', 'SMARCA4'), ('TP53', 'SMARCB1'), ('TP53', 'SMARCB1'), ('TP53', 'SMARCB1'), ('TP53', 'SMARCA4'), ('TP53', 'SMARCC1'), ('TP53', 'UBE2I'), ('TP53', 'CHD3'), ('RUVBL1', 'TRRAP'), ('RUVBL1', 'SMARCA2'), ('RUVBL1', 'ACTL6A'), ('RUVBL1', 'ACTL6A'), ('SRCAP', 'NCOR1'), ('TP53', 'MIF'), ('TP53', 'S100A8'), ('TP53', 'CCT5'), ('TP53', 'HSPB1'), ('TP53', 'WDR33'), ('TP53', 'EIF2S2'), ('TP53', 'SNRPN'), ('TP53', 'ERH'), ('TP53', 'HSP90AA1'), ('POLR2E', 'POLR2A'), ('POLR2E', 'POLR2B'), ('POLR2E', 'POLR2C'), ('POLR2E', 'POLR2A'), ('POLR2E', 'POLR2B'), ('POLR2E', 'POLR2C'), ('POLR2E', 'POLR2E'), ('POLR2E', 'POLR2H'), ('POLR2E', 'POLR2H'), ('POLR2E', 'TAF15'), ('POLR2E', 'TAF15'), ('TP53', 'TOP1'), ('TP53', 'PRKDC'), ('TP53', 'HMGB1'), ('TP53', 'HMGB1'), ('TP53', 'HMGB1'), ('TP53', 'HMGB1'), ('TP53', 'APEX1'), ('TP53', 'RFC1'), ('TP53', 'CDC2'), ('TP53', 'RPA1'), ('TP53', 'YWHAZ'), ('TP53', 'YWHAZ'), ('TP53', 'CDC2'), ('TP53', 'ZNF148'), ('TP53', 'ZNF148'), ('TP53', 'MED1'), ('POLR2E', 'XRCC5'), ('TP53', 'MED1'), ('TP53', 'UBE2I'), ('TP53', 'TOP1'), ('TP53', 'CEBPZ'), ('TP53', 'CEBPZ'), ('TP53', 'HMGB1'), ('TP53', 'NCL'), ('TP53', 'NCL'), ('TP53', 'NCL'), ('TP53', 'PARP1'), ('TP53', 'UBE2I'), ('RUVBL2', 'PCNA'), ('TP53', 'TOP2A'), ('TP53', 'TOP2B'), ('TP53', 'TOP2A'), ('TP53', 'TOP2B'), ('TP53', 'GNL3'), ('TP53', 'GNL3'), ('TP53', 'GNL3'), ('TP53', 'USP7'), ('TP53', 'USP7'), ('TP53', 'USP7'), ('TP53', 'USP7'), ('TP53', 'HSPA9'), ('TP53', 'HSPA9'), ('TP53', 'TP53BP1'), ('TP53', 'TP53BP1'), ('TP53', 'TP53BP1'), ('TP53', 'KPNB1'), ('TP53', 'KPNB1'), ('TP53', 'HSP90AA1'), ('TP53', 'HSPA8'), ('TP53', 'UBE2I'), ('NUFIP1', 'FMR1'), ('NUFIP1', 'FMR1'), ('TP53', 'YBX1'), ('TP53', 'YBX1'), ('TP53', 'HSP90AA1'), ('TP53', 'FKBP3'), ('TP53', 'CCAR1'), ('TP53', 'HSPA4'), ('TP53', 'STUB1'), ('TP53', 'SERPINH1'), ('TP53', 'VRK1'), ('RUVBL1', 'ING3'), ('RUVBL2', 'ING3'), ('RUVBL1', 'ING3'), ('TP53', 'SFN'), ('TP53', 'TXN'), ('TP53', 'USP7'), ('TP53', 'USP7'), ('TP53', 'USP39'), ('TFPT', 'UCHL5'), ('INO80E', 'UCHL5'), ('RUVBL1', 'MEPCE'), ('RUVBL2', 'MEPCE'), ('PIH1D1', 'RPAP2'), ('RUVBL2', 'RPAP2'), ('PIH1D1', 'RPAP3'), ('RUVBL1', 'PCBP2'), ('RUVBL1', 'RUVBL2'), ('RUVBL2', 'RUVBL1'), ('RUVBL2', 'DDX42'), ('RUVBL2', 'RPAP3'), ('RUVBL2', 'KIAA1967'), ('RUVBL2', 'MATR3'), ('PIH1D1', 'RUVBL2'), ('RUVBL2', 'PIH1D1'), ('RUVBL2', 'PCBP2'), ('RUVBL2', 'POLR2B'), ('RUVBL2', 'RFC4'), ('RP11-529I10.4', 'RUVBL2'), ('RUVBL2', 'RP11-529I10.4'), ('RUVBL1', 'RUVBL2'), ('RUVBL2', 'RUVBL1'), ('SRCAP', 'RUVBL2'), ('RUVBL2', 'SRCAP'), ('UXT', 'RUVBL2'), ('RUVBL2', 'UXT'), ('RUVBL2', 'ACTB'), ('RUVBL2', 'C19orf2'), ('POLR2E', 'RPAP3'), ('RUVBL2', 'YEATS4'), ('RUVBL2', 'EFTUD2'), ('RUVBL2', 'C1orf57'), ('ZNHIT2', 'RUVBL2'), ('RUVBL2', 'ZNHIT2'), ('RUVBL2', 'PFDN6'), ('RUVBL1', 'HSPA8'), ('RUVBL2', 'PDRG1'), ('RUVBL1', 'C1orf57'), ('POLR2E', 'RUVBL2'), ('RUVBL2', 'POLR2E'), ('RUVBL1', 'ACTB'), ('RUVBL2', 'KPNB1'), ('RUVBL2', 'HSPA5'), ('RUVBL2', 'PPP1CA'), ('RUVBL2', 'HNRNPU'), ('RUVBL1', 'DDX42'), ('RUVBL2', 'RFC5'), ('RUVBL2', 'RPS4X'), ('RUVBL2', 'POLR2H'), ('RP11-529I10.4', 'RUVBL1'), ('RUVBL1', 'RP11-529I10.4'), ('RUVBL2', 'PCBP1'), ('RUVBL2', 'ACTL6A'), ('UXT', 'RPAP3'), ('RUVBL2', 'TRIM28'), ('RUVBL1', 'KPNB1'), ('RUVBL2', 'DNM2'), ('UXT', 'RPAP2'), ('PIH1D1', 'RUVBL1'), ('RUVBL1', 'PIH1D1'), ('ACTR6', 'RUVBL2'), ('RUVBL2', 'ACTR6'), ('RUVBL1', 'HSPA1L'), ('RUVBL2', 'NOP58'), ('RUVBL2', 'NFRKB'), ('TFPT', 'RUVBL2'), ('RUVBL2', 'TFPT'), ('RUVBL1', 'ELAVL1'), ('RUVBL2', 'DMAP1'), ('INO80E', 'RUVBL2'), ('RUVBL2', 'INO80E'), ('RUVBL1', 'DMAP1'), ('RUVBL2', 'UCHL5'), ('RUVBL1', 'YEATS4'), ('RUVBL2', 'RPL38'), ('RUVBL2', 'ACTR5'), ('ACTR5', 'RUVBL2'), ('RUVBL1', 'PCBP1'), ('POLR2E', 'TAF10'), ('RUVBL2', 'RPS2'), ('ACTR8', 'RUVBL2'), ('RUVBL2', 'ACTR8'), ('RUVBL2', 'HSPA1L'), ('TP53', 'USP7'), ('TP53', 'UBE2N'), ('TP53', 'UBE2N'), ('TP53', 'UBE2N'), ('TP53', 'HDAC2'), ('C20orf20', 'BRD8'), ('C20orf20', 'BRD8'), ('C20orf20', 'BRD8'), ('TP53', 'BRD8'), ('TP53', 'DNAJB6'), ('TP53', 'EIF2C2'), ('TP53', 'IGF2BP1'), ('TP53', 'LRPPRC'), ('TP53', 'SMARCD2'), ('RUVBL1', 'CTNNB1')]
data1 = open('/Users/petercimermancic/Documents/Krogan_lab/FalsePositivies.txt')
DB_FP = [d.strip().split('\t')[2] for d in data1.readlines()]
data1.close()


#saint,nsaf,zscore,dscore,mist = {},{},{},{},{}

from operator import itemgetter

SAINT,NSAF,ZSCORE,DSCORE,MIST = [],[],[],[],[]
XSAINT,XNSAF,XZSCORE,XDSCORE,XMIST = [],[],[],[],[]

for x in xrange(300,len(mist),500):
    
    S = sorted(saint.values(),reverse=True)[x]
    SAINT.append(len( set([s for s in saint if saint[s] >= S]) & set(DB) ) / float(x))
    XSAINT.append(len([s for s in saint if saint[s] >= S]))

    N = sorted(nsaf.values(),reverse=True)[x]
    NSAF.append(len( set([s for s in nsaf if nsaf[s] >= N]) & set(DB) ) / float(x))
    XNSAF.append(len([s for s in nsaf if nsaf[s] >= N]))

    Z = sorted(zscore.values(),reverse=True)[x]
    ZSCORE.append(len( set([s for s in zscore if zscore[s] >= Z]) & set(DB) ) / float(x))
    XZSCORE.append(len([s for s in zscore if zscore[s] >= Z]))

    D = sorted(dscore.values(),reverse=True)[x]
    DSCORE.append(len( set([s for s in dscore if dscore[s] >= D]) & set(DB) ) / float(x))
    XDSCORE.append(len([s for s in dscore if dscore[s] >= D]))

    M = sorted(mist.values(),reverse=True)[x]
    MIST.append(len( set([s for s in mist if mist[s] >= M]) & set(DB) ) / float(x))
    XMIST.append(len([s for s in mist if mist[s] >= M]))

X = range(300,len(mist),300)
print len(DB),len(saint)


import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_title('TIP49 dataset')
ax.plot( XSAINT,SAINT,'r--', XNSAF,NSAF,'c--', XZSCORE,ZSCORE,'m--', XDSCORE,DSCORE,'g--', XMIST,MIST,'b--',)
ax.legend(('SAInt','PP-NSAF','Z-score','Weighted D-score','MiST'))

ax.plot( XSAINT,SAINT,'ro', XZSCORE,NSAF,'co', XZSCORE,ZSCORE,'mo', XDSCORE,DSCORE,'go', XMIST,MIST,'bo',)

#plt.savefig('DatabaseComparison.svg',format='SVG')
plt.show()
