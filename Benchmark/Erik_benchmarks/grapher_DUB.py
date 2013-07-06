import pickle
from operator import itemgetter



DB = [('USP11', 'RANBP9'), ('USP16', 'EXOSC10'), ('USP53', 'TRAF2'), ('USP2', 'HOOK2'), ('USP2', 'C1orf165'), ('STAMBPL1', 'TRIP13'), ('USP2', 'TRAF2'), ('EIF3S5', 'HAX1'), ('USP5', 'TADA3L'), ('USP5', 'TADA3L'), ('USP11', 'TCEAL1'), ('USP8', 'YWHAB'), ('EIF3S3', 'EIF4A2'), ('PSMD14', 'WDR71'), ('PSMD7', 'WDR71'), ('EIF3S3', 'EIF3S10'), ('PSMD7', 'PSMD10'), ('PSMD14', 'PSMD13'), ('PSMD7', 'PSMD13'), ('EIF3S5', 'EIF3S10'), ('USP39', 'CD2BP2'), ('USP8', 'YWHAZ'), ('PSMD14', 'PSMD6'), ('PSMD7', 'PSMD6'), ('USP7', 'RUVBL2'), ('USP15', 'PHB2'), ('BRCC3', 'BRE'), ('USP8', 'YWHAQ'), ('PSMD7', 'PSMC4'), ('PARP11', 'NUP98'), ('USP7', 'RNPS1'), ('BRCC3', 'BRE'), ('BRCC3', 'BRE'), ('BRCC3', 'FANCB'), ('BRCC3', 'FANCB'), ('BRCC3', 'FANCB'), ('BRCC3', 'TP53'), ('BRCC3', 'TP53'), ('BAP1', 'SFN'), ('USP8', 'SFN'), ('BRCC3', 'BRE'), ('BRCC3', 'TP53'), ('USP8', 'OTUB1'), ('USP52', 'PAN3'), ('USP52', 'PABPC1'), ('EIF3S5', 'EIF3S10'), ('EIF3S3', 'EIF3S10'), ('ATXN3', 'VCP'), ('COPS5', 'S100A7'), ('COPS5', 'S100A7'), ('USP33', 'DDX21'), ('USP33', 'DDX56'), ('COPS6', 'FAU'), ('USP33', 'ITGB4BP'), ('COPS6', 'RPL15'), ('UCHL1', 'CBX1'), ('EIF3S5', 'DKC1'), ('EIF3S5', 'PCID1'), ('COPS6', 'EMD'), ('COPS6', 'ERH'), ('COPS6', 'HMOX2'), ('COPS6', 'MYCBP'), ('UCHL1', 'NEDD8'), ('COPS6', 'PFKL'), ('COPS6', 'SMN1'), ('USP13', 'MYO15A'), ('COPS6', 'BTBD2'), ('COPS6', 'MAP7D1'), ('COPS6', 'PRKRA'), ('COPS6', 'PSMD11'), ('COPS6', 'RFC5'), ('COPS6', 'RPA2'), ('COPS6', 'SNRPG'), ('COPS6', 'TK1'), ('COPS6', 'ZNF24'), ('EIF3S5', 'EEF1A1'), ('UCHL1', 'RANBP9'), ('COPS5', 'UCHL1'), ('UCHL1', 'COPS5'), ('COPS5', 'UCHL1'), ('UCHL1', 'COPS5'), ('COPS5', 'COPS8'), ('COPS6', 'COPS8'), ('BAP1', 'UBE2D2'), ('BAP1', 'WBP2'), ('USP8', 'GRB2'), ('USP10', 'G3BP1'), ('USP10', 'G3BP1'), ('TNFAIP3', 'TRAF2'), ('TNFAIP3', 'TRAF2'), ('TNFAIP3', 'TNFAIP3'), ('EIF3S5', 'FRAP1'), ('EIF3S5', 'FRAP1'), ('EIF3S5', 'FRAP1'), ('UCHL3', 'NEDD8'), ('UCHL3', 'NEDD8'), ('COPS5', 'MIF'), ('COPS5', 'MIF'), ('TNFAIP3', 'IKBKG'), ('TNFAIP3', 'IKBKG'), ('BAP1', 'NET1'), ('COPS6', 'EIF3S6'), ('COPS6', 'EIF3S6'), ('TNFAIP3', 'YWHAZ'), ('TNFAIP3', 'YWHAH'), ('TNFAIP3', 'YWHAH'), ('BAP1', 'ACTN4'), ('BAP1', 'ACTN4'), ('USP11', 'RANBP9'), ('USP11', 'RANBP9'), ('ATXN3', 'RAD23A'), ('ATXN3', 'RAD23B'), ('ATXN3', 'RAD23A'), ('ATXN3', 'RAD23B'), ('ATXN3', 'RAD23A'), ('TNFAIP3', 'YWHAB'), ('TNFAIP3', 'YWHAH'), ('TNFAIP3', 'YWHAZ'), ('TNFAIP3', 'YWHAB'), ('TNFAIP3', 'YWHAH'), ('TNFAIP3', 'YWHAZ'), ('TNFAIP3', 'YWHAB'), ('TNFAIP3', 'YWHAH'), ('TNFAIP3', 'YWHAZ'), ('TNFAIP3', 'YWHAE'), ('TNFAIP3', 'IKBKB'), ('TNFAIP3', 'CHUK'), ('USP7', 'TP53'), ('USP7', 'TP53'), ('USP7', 'TP53'), ('USP7', 'TP53'), ('USP33', 'TCEB1'), ('USP7', 'TRAF2'), ('USP7', 'TRAF3'), ('EIF3S5', 'EIF3S10'), ('EIF3S3', 'EIF3S10'), ('EIF3S3', 'EIF3S4'), ('BAP1', 'BAI1'), ('BAP1', 'BAI1'), ('PSMD7', 'ATXN7'), ('USP33', 'SELENBP1'), ('USP33', 'SELENBP1'), ('OTUB1', 'DDX5'), ('OTUB1', 'GNB2L1'), ('OTUB1', 'NPM1'), ('OTUB1', 'EBNA1BP2'), ('OTUB1', 'EIF4A3'), ('OTUB1', 'DDX54'), ('OTUB1', 'DDX23'), ('OTUB1', 'PCNA'), ('OTUB1', 'NAT10'), ('OTUB1', 'DDX24'), ('OTUB1', 'GNB2L1'), ('OTUB1', 'FUS'), ('USP4', 'TRIM21'), ('UCHL1', 'TP53'), ('CYLD', 'TBK1'), ('USP8', 'KIF23'), ('USP22', 'ATXN7L3'), ('USP22', 'TRRAP'), ('USP22', 'ATXN7'), ('USP22', 'GCN5L2'), ('USP22', 'ATXN7L3'), ('USP22', 'SUPT3H'), ('USP22', 'TAF10'), ('USP22', 'ATXN7'), ('USP22', 'ENY2'), ('ATXN3', 'VCP'), ('COPS5', 'CPSF1'), ('OTUD5', 'TRAF3'), ('USP7', 'PDCD6IP'), ('COPS5', 'TP53'), ('COPS5', 'TP53'), ('COPS5', 'TP53'), ('PSMD14', 'PSMD10'), ('PSMD14', 'PSMA3'), ('PSMD14', 'PSMC6'), ('PSMD14', 'PSMD4'), ('USP7', 'GMPS'), ('USP7', 'GMPS'), ('USP7', 'GMPS'), ('USP7', 'TP53'), ('USP7', 'TP53'), ('BRCC3', 'UIMC1'), ('CYLD', 'PLK1'), ('CYLD', 'PLK1'), ('USP39', 'BCDIN3'), ('UCHL5', 'RUVBL2'), ('COPS5', 'COPS7A'), ('USP7', 'TP53'), ('STAMBP', 'VPS24'), ('STAMBP', 'VPS24'), ('STAMBP', 'VPS24'), ('UCHL5', 'ADRM1'), ('UCHL5', 'ADRM1'), ('UCHL5', 'PSMD2'), ('UCHL5', 'ADRM1'), ('UCHL5', 'ADRM1'), ('COPS5', 'DDB1'), ('COPS6', 'DDB1'), ('COPS5', 'CUL4A'), ('COPS6', 'CUL4A'), ('COPS5', 'DDB1'), ('COPS6', 'DDB1'), ('COPS5', 'TRAF2'), ('COPS5', 'TRAF2'), ('COPS5', 'TP53'), ('UCHL5', 'ADRM1'), ('UCHL5', 'ADRM1'), ('UCHL5', 'ADRM1'), ('UCHL5', 'ADRM1'), ('USP28', 'TP53BP1'), ('USP28', 'TP53BP1'), ('USP28', 'TP53BP1'), ('STAMBP', 'VPS24'), ('STAMBP', 'VPS24'), ('BAP1', 'UBE2D2'), ('BAP1', 'UBE2V1'), ('USP25', 'GAPDH'), ('USP21', 'GAPDH'), ('USP25', 'ACTA1'), ('USP21', 'ACTA1'), ('USP25', 'FLNC'), ('USP21', 'FLNC'), ('USP25', 'ACTA1'), ('USP21', 'ACTA1'), ('USP4', 'TRIM21'), ('STAMBP', 'CLTC'), ('STAMBPL1', 'CLTC'), ('STAMBP', 'CLTC'), ('STAMBPL1', 'CLTC')]

data1 = open('DUB.txt','r')
D1 = data1.readlines()
data1.close()

comppass = []
saint = []
for d in D1[1:]:
    d = d.strip().split('\t')
    comppass.append( (d[0],d[1],float(d[6])) )
    saint.append( (d[0],d[1],float(d[7])) )

data2 = open('Dub_mist.txt','r')
D2 = data2.readlines()
data2.close()

mist = []
for d in D2:
    d = d.strip().split('\t')
    mist.append( (d[0],d[1],float(d[2])) )

SAINT,COMPPASS,MIST = [],[],[]
XSAINT,XCOMPPASS,XMIST = [],[],[]

for x in xrange(800,10000,1000):
    
    S = sorted(saint,key=itemgetter(2),reverse=True)[x][2]
    SAINT.append(len( set([(s[0],s[1]) for s in saint if s[2] >= S]) & set(DB) ) / float(x))
    XSAINT.append(len([(s[0],s[1]) for s in saint if s[2] >= S]))

    N = sorted(comppass,key=itemgetter(2),reverse=True)[x][2]
    COMPPASS.append(len( set([(s[0],s[1]) for s in comppass if s[2] >= N]) & set(DB) ) / float(x))
    XCOMPPASS.append(len([(s[0],s[1]) for s in comppass if s[2] >= N]))

    M = sorted(mist,key=itemgetter(2),reverse=True)[x][2]
    MIST.append(len( set([(s[0],s[1]) for s in mist if s[2] >= M]) & set(DB) ) / float(x))
    XMIST.append(len([(s[0],s[1]) for s in mist if s[2] >= M]))

X = range(300,len(mist),300)
print len(DB),len(saint)

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_title('DUB dataset')
ax.plot( XSAINT[:-2],SAINT[:-2],'r--', XCOMPPASS,COMPPASS,'c--', XMIST,MIST,'b--',)
ax.legend(('SAInt','CompPASS','MiST'))

ax.plot( XSAINT[:-2],SAINT[:-2],'ro', XCOMPPASS,COMPPASS,'co', XMIST,MIST,'bo',)

plt.savefig('DatabaseComparison2.svg',format='SVG')
#plt.show()

