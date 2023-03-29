import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#--------------------------------------------------------------------
# begin user input
#--------------------------------------------------------------------

solute_list = ['Na', 'K', 'Cl', 'HCO3', 'NH4', 'Volume']

segment_early = ['PT', 'S3', 'SDL', 'mTAL', 'cTAL', 'DCT', 'CNT']
segment_jux = ['SDL', 'LDL', 'LAL']
segment_cd = ['CCD','OMCD','IMCD']

compare = 3 #2, 3, or 4

direct1 = 'NP-baseline2022-04-22'
sex1 = 'female'

direct2 = 'midpregnant_rat'
sex2 = 'midpregnant'

direct3 = 'latepregnant_rat'
sex3 = 'latepregnant'

direct4 = ''
sex4 = ''

label1 = 'virgin'
label2 = 'MP'
label3 = 'LP'
label4 = ''

save_tiff = 0 # set to 1 if want to save .tiff file
humOrrat = 'rat'
# conversion factors
if humOrrat == 'rat':
    sup_ratio = 2.0/3.0
    neph_per_kidney = 36000 #number of nephrons per kidney
elif humOrrat == 'hum':
    sup_ratio = 0.85
    neph_per_kidney = 1000000
jux_ratio = 1-sup_ratio
neph_weight = [sup_ratio, 0.4*jux_ratio, 0.3*jux_ratio, 0.15*jux_ratio, 0.1*jux_ratio, 0.05*jux_ratio ]


p_to_mu = 1e-6 #convert pmol to micromole, same for nl to ml for volume
cf = neph_per_kidney * p_to_mu

mu_conv = 1e-3 #from mumol to mmol
min_per_day = 1440
sol_conv = mu_conv * min_per_day
vol_conv = min_per_day

segment_transport = ['PT','DL','LAL','TAL','DCT','CNT','CD']


#-----------------------------------------------------------------------
#end user input
#----------------------------------------------------------------------
#========================================
# functions used
#========================================

def get_transport(direct, sex, solute, segment, supOrjux):
    os.chdir(direct)
    if solute == 'Volume':
        fname = sex+'_'+humOrrat+'_'+segment+'_water_volume_in_Lumen'+supOrjux+'.txt'
        file = open(fname, 'r')
        start = float(file.readline())
        end = float(np.loadtxt(fname, delimiter = '\n', unpack = True)[-1])
        file.close()
        transport = start - end
    else:
        fname = sex+'_'+humOrrat+'_'+segment+'_flow_of_'+solute+'_in_Lumen'+supOrjux+'.txt'
        file = open(fname, 'r')
        start = float(file.readline())
        end = float(np.loadtxt(fname, delimiter = '\n', unpack = True)[-1])
        file.close()
        transport = start - end
    os.chdir('..')
    transport_cf = transport * cf
    return transport_cf

def get_trans_data(direct, sex, solute, segments):
    # get transport data for nephron segments (i.e., not collecting duct)
    sup_trans = np.zeros(len(segments))
    jux1_trans = np.zeros(len(segments))
    jux2_trans = np.zeros(len(segments))
    jux3_trans = np.zeros(len(segments))
    jux4_trans = np.zeros(len(segments))
    jux5_trans = np.zeros(len(segments))
    
    for s in range(len(segments)):
        seg = segments[s]
        if seg[-2:].lower() == 'cd':
            print('segment: ' + seg)
            raise Exception('not for collecting duct, use get_cd_data for cd')
        if seg.lower() != 'ldl' and seg.lower() != 'lal':
            sup_trans[s] = get_transport(direct, sex, solute, seg, '_sup')
        jux1_trans[s] = get_transport(direct, sex, solute, seg, '_jux1')
        jux2_trans[s] = get_transport(direct, sex, solute, seg, '_jux2')
        jux3_trans[s] = get_transport(direct, sex, solute, seg, '_jux3')
        jux4_trans[s] = get_transport(direct, sex, solute, seg, '_jux4')
        jux5_trans[s] = get_transport(direct, sex, solute, seg, '_jux5')
        
    trans_vals_weighted = np.matrix([sup_trans*neph_weight[0], jux1_trans*neph_weight[1],
                                     jux2_trans*neph_weight[2], jux3_trans*neph_weight[3],
                                     jux4_trans*neph_weight[4], jux5_trans*neph_weight[5]])
        
    sup_vals = neph_weight[0]*sup_trans
    jux_vals = neph_weight[1]*jux1_trans + neph_weight[2]*jux2_trans + \
        neph_weight[3]*jux3_trans + neph_weight[4]*jux4_trans + neph_weight[5]*jux5_trans
    
    return sup_vals, jux_vals, trans_vals_weighted

def get_cd_data(direct, sex, solute, segments):
    trans = np.zeros(len(segments))
    
    for s in range(len(segments)):
        seg = segments[s]
        if seg[-2:].lower() != 'cd':
            print('segment: ' + seg)
            raise Exception('only for collecting duct segments')
        trans[s] = get_transport(direct, sex, solute, seg, '')
    trans = trans
    return trans

def get_data2(solute):
    #segment_transport = ['PT','DL','LAL','TAL','DCT','CNT','CD']
    sup1, jux1, trans1 = get_trans_data(direct1, sex1, solute, segment_early)
    sup2, jux2, trans2 = get_trans_data(direct2, sex2, solute, segment_early)
    
    dl1 = get_trans_data(direct1, sex1, solute, segment_jux)[1]
    dl2 = get_trans_data(direct2, sex2, solute, segment_jux)[1]
    
    cd1 = get_cd_data(direct1, sex1, solute, segment_cd)
    cd2 = get_cd_data(direct2, sex2, solute, segment_cd)
    
    sup1_final = np.zeros(len(segment_transport))
    jux1_final = np.zeros(len(segment_transport))
    sup2_final = np.zeros(len(segment_transport))
    jux2_final = np.zeros(len(segment_transport))
    
    # PT = pct & s3
    sup1_final[0] = sum(sup1[0:1+1]) 
    sup2_final[0] = sum(sup2[0:1+1])
    
    jux1_final[0] = sum(jux1[0:1+1])
    jux2_final[0] = sum(jux2[0:1+1])
    
    # DL = sdl for sup, sdl + ldl for jux
    sup1_final[1] = sup1[2]
    sup2_final[1] = sup2[2]
    
    #sdl + ldl
    jux1_final[1] = sum(dl1[0:1+1])
    jux2_final[1] = sum(dl2[0:1+1])
    
    # LAL
    jux1_final[2] = dl1[-1]
    jux2_final[2] = dl2[-1]
    
    # TAL = mTAL + cTAL
    sup1_final[3] = sum(sup1[3:4+1])
    sup2_final[3] = sum(sup2[3:4+1])
    
    jux1_final[3] = sum(jux1[3:4+1])
    jux2_final[3] = sum(jux2[3:4+1])
    
    # DCT, CNT
    sup1_final[4:5+1] = sup1[5:6+1]
    sup2_final[4:5+1] = sup2[5:6+1]
    
    jux1_final[4:5+1] = jux1[5:6+1]
    jux2_final[4:5+1] = jux2[5:6+1]
    
    # CD = ccd + omcd + imcd
    sup1_final[-1] = sum(cd1)
    sup2_final[-1] = sum(cd2)
    
    vals1 = np.zeros(2*7).reshape(2,7)
    vals1[0] = sup1_final
    vals1[1] = jux1_final
    
    vals2 = np.zeros(2*7).reshape(2,7)
    vals2[0] = sup2_final
    vals2[1] = jux2_final
    
    
    return vals1, vals2

def get_data3(solute):
    #segment_transport = ['PT','DL','LAL','TAL','DCT','CNT','CD']
    sup1, jux1, trans1 = get_trans_data(direct1, sex1, solute, segment_early)
    sup2, jux2, trans2 = get_trans_data(direct2, sex2, solute, segment_early)
    sup3, jux3, trans3 = get_trans_data(direct3, sex3, solute, segment_early)
    
    dl1 = get_trans_data(direct1, sex1, solute, segment_jux)[1]
    dl2 = get_trans_data(direct2, sex2, solute, segment_jux)[1]
    dl3 = get_trans_data(direct3, sex3, solute, segment_jux)[1]
    
    cd1 = get_cd_data(direct1, sex1, solute, segment_cd)
    cd2 = get_cd_data(direct2, sex2, solute, segment_cd)
    cd3 = get_cd_data(direct3, sex3, solute, segment_cd)
    
    sup1_final = np.zeros(len(segment_transport))
    jux1_final = np.zeros(len(segment_transport))
    sup2_final = np.zeros(len(segment_transport))
    jux2_final = np.zeros(len(segment_transport))
    sup3_final = np.zeros(len(segment_transport))
    jux3_final = np.zeros(len(segment_transport))
    
    # PT = pct & s3
    sup1_final[0] = sum(sup1[0:1+1]) 
    sup2_final[0] = sum(sup2[0:1+1])
    sup3_final[0] = sum(sup3[0:1+1])
    
    jux1_final[0] = sum(jux1[0:1+1])
    jux2_final[0] = sum(jux2[0:1+1])
    jux3_final[0] = sum(jux3[0:1+1])
    
    # DL = sdl for sup, sdl + ldl for jux
    sup1_final[1] = sup1[2]
    sup2_final[1] = sup2[2]
    sup3_final[1] = sup3[2]
    
    #sdl + ldl
    jux1_final[1] = sum(dl1[0:1+1])
    jux2_final[1] = sum(dl2[0:1+1])
    jux3_final[1] = sum(dl3[0:1+1])
    
    # LAL
    jux1_final[2] = dl1[-1]
    jux2_final[2] = dl2[-1]
    jux3_final[2] = dl3[-1]
    
    # TAL = mTAL + cTAL
    sup1_final[3] = sum(sup1[3:4+1])
    sup2_final[3] = sum(sup2[3:4+1])
    sup3_final[3] = sum(sup3[3:4+1])
    
    jux1_final[3] = sum(jux1[3:4+1])
    jux2_final[3] = sum(jux2[3:4+1])
    jux3_final[3] = sum(jux3[3:4+1])
    
    # DCT, CNT
    sup1_final[4:5+1] = sup1[5:6+1]
    sup2_final[4:5+1] = sup2[5:6+1]
    sup3_final[4:5+1] = sup3[5:6+1]
    
    jux1_final[4:5+1] = jux1[5:6+1]
    jux2_final[4:5+1] = jux2[5:6+1]
    jux3_final[4:5+1] = jux3[5:6+1]
    
    # CD = ccd + omcd + imcd
    sup1_final[-1] = sum(cd1)
    sup2_final[-1] = sum(cd2)
    sup3_final[-1] = sum(cd3)
    
    vals1 = np.zeros(2*7).reshape(2,7)
    vals1[0] = sup1_final
    vals1[1] = jux1_final
    
    vals2 = np.zeros(2*7).reshape(2,7)
    vals2[0] = sup2_final
    vals2[1] = jux2_final
    
    vals3 = np.zeros(2*7).reshape(2,7)
    vals3[0] = sup3_final
    vals3[1] = jux3_final
    
    
    return vals1, vals2, vals3


def get_data4(solute):
    #segment_transport = ['PT','DL','LAL','TAL','DCT','CNT','CD']
    sup1, jux1, trans1 = get_trans_data(direct1, sex1, solute, segment_early)
    sup2, jux2, trans2 = get_trans_data(direct2, sex2, solute, segment_early)
    sup3, jux3, trans3 = get_trans_data(direct3, sex3, solute, segment_early)
    sup4, jux4, trans4 = get_trans_data(direct4, sex4, solute, segment_early)
    
    dl1 = get_trans_data(direct1, sex1, solute, segment_jux)[1]
    dl2 = get_trans_data(direct2, sex2, solute, segment_jux)[1]
    dl3 = get_trans_data(direct3, sex3, solute, segment_jux)[1]
    dl4 = get_trans_data(direct4, sex4, solute, segment_jux)[1]
    
    cd1 = get_cd_data(direct1, sex1, solute, segment_cd)
    cd2 = get_cd_data(direct2, sex2, solute, segment_cd)
    cd3 = get_cd_data(direct3, sex3, solute, segment_cd)
    cd4 = get_cd_data(direct4, sex4, solute, segment_cd)
    
    sup1_final = np.zeros(len(segment_transport))
    jux1_final = np.zeros(len(segment_transport))
    sup2_final = np.zeros(len(segment_transport))
    jux2_final = np.zeros(len(segment_transport))
    sup3_final = np.zeros(len(segment_transport))
    jux3_final = np.zeros(len(segment_transport))
    sup4_final = np.zeros(len(segment_transport))
    jux4_final = np.zeros(len(segment_transport))
    
    # PT = pct & s3
    sup1_final[0] = sum(sup1[0:1+1]) 
    sup2_final[0] = sum(sup2[0:1+1])
    sup3_final[0] = sum(sup3[0:1+1])
    sup4_final[0] = sum(sup4[0:1+1])
    
    jux1_final[0] = sum(jux1[0:1+1])
    jux2_final[0] = sum(jux2[0:1+1])
    jux3_final[0] = sum(jux3[0:1+1])
    jux4_final[0] = sum(jux4[0:1+1])
    
    # DL = sdl for sup, sdl + ldl for jux
    sup1_final[1] = sup1[2]
    sup2_final[1] = sup2[2]
    sup3_final[1] = sup3[2]
    sup4_final[1] = sup4[2]
    
    #sdl + ldl
    jux1_final[1] = sum(dl1[0:1+1])
    jux2_final[1] = sum(dl2[0:1+1])
    jux3_final[1] = sum(dl3[0:1+1])
    jux4_final[1] = sum(dl4[0:1+1])
    
    # LAL
    jux1_final[2] = dl1[-1]
    jux2_final[2] = dl2[-1]
    jux3_final[2] = dl3[-1]
    jux4_final[2] = dl4[-1]
    
    # TAL = mTAL + cTAL
    sup1_final[3] = sum(sup1[3:4+1])
    sup2_final[3] = sum(sup2[3:4+1])
    sup3_final[3] = sum(sup3[3:4+1])
    sup4_final[3] = sum(sup4[3:4+1])
    
    jux1_final[3] = sum(jux1[3:4+1])
    jux2_final[3] = sum(jux2[3:4+1])
    jux3_final[3] = sum(jux3[3:4+1])
    jux4_final[3] = sum(jux4[3:4+1])
    
    # DCT, CNT
    sup1_final[4:5+1] = sup1[5:6+1]
    sup2_final[4:5+1] = sup2[5:6+1]
    sup3_final[4:5+1] = sup3[5:6+1]
    sup4_final[4:5+1] = sup4[5:6+1]
    
    jux1_final[4:5+1] = jux1[5:6+1]
    jux2_final[4:5+1] = jux2[5:6+1]
    jux3_final[4:5+1] = jux3[5:6+1]
    jux4_final[4:5+1] = jux4[5:6+1]
    
    # CD = ccd + omcd + imcd
    sup1_final[-1] = sum(cd1)
    sup2_final[-1] = sum(cd2)
    sup3_final[-1] = sum(cd3)
    sup4_final[-1] = sum(cd4)
    
    vals1 = np.zeros(2*7).reshape(2,7)
    vals1[0] = sup1_final
    vals1[1] = jux1_final
    
    vals2 = np.zeros(2*7).reshape(2,7)
    vals2[0] = sup2_final
    vals2[1] = jux2_final
    
    vals3 = np.zeros(2*7).reshape(2,7)
    vals3[0] = sup3_final
    vals3[1] = jux3_final
    
    vals4 = np.zeros(2*7).reshape(2,7)
    vals4[0] = sup4_final
    vals4[1] = jux4_final
    
    return vals1, vals2, vals3, vals4
    
                    
#===============================================================================
#   make the figures
#===============================================================================
# colors
c1 = 'c'
c2 = 'mediumvioletred'
c3 = 'green'
c4 = 'purple'
juxc = 'white'
al = 0.4

title_size = 25
ylab_size = 20
leg_size=18
fig_lab_size = 25
xlab_size=20
in_xlab_size=16

# figure specs
figlab_shift = -1.5
width = 0.25


fig = plt.figure(figsize = (20,20))
plt.rcParams.update({'font.size':20})


# positiions
x = np.arange(len(segment_transport))


# positiions
sup_pos = np.arange(len(segment_transport[:6]))
full_pos = np.arange(len(segment_transport))
later_pos = np.arange(len(segment_transport[:6]), len(segment_transport))


#------------
# Na
#------------
if compare == 2:
    vals1, vals2 = get_data2('Na')
elif compare == 3:
    vals1, vals2, vals3 = get_data3('Na')
else:
    vals1, vals2, vals3, vals4 = get_data4('Na')


ax1 = fig.add_subplot(3,2,1)

# bar1
sup1 = ax1.bar(full_pos-width, vals1[0], width, align = 'center', edgecolor = 'k', color = c1, label = label1)
jux1 = ax1.bar(full_pos-width, vals1[1], width, bottom = vals1[0], align = 'center', edgecolor = 'k', color = c1, alpha = al)
# bar2
sup2 = ax1.bar(full_pos, vals2[0], width, align = 'center', edgecolor = 'k', color = c2, label = label2)
jux2 = ax1.bar(full_pos, vals2[1], width, bottom = vals2[0], align = 'center', edgecolor = 'k', color = c2, alpha = al)

if compare>2:
    # bar 3
    sup3 = ax1.bar(full_pos+width, vals3[0], width, align = 'center', edgecolor = 'k', color = c3, label = label3)
    jux3 = ax1.bar(full_pos+width, vals3[1], width, bottom = vals3[0], align = 'center', edgecolor = 'k', color = c3, alpha = al)
    if compare >3:
        # bar 4
        sup4 = ax1.bar(full_pos+ 2*width, vals4[0], width, align = 'center', edgecolor = 'k', color = c4, label = label4)
        jux4 = ax1.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], align = 'center', edgecolor = 'k', color = c4, alpha = al)
plt.axhline(0, color = 'k')
# inset
axins1 = inset_axes(ax1, width="35%", height = 1.0, loc = 5)
axins1.bar(full_pos-width, vals1[0], width, color = c1, edgecolor = 'k')
axins1.bar(full_pos-width, vals1[1], width, bottom = vals1[0], color = c1, alpha = al, edgecolor = 'k')
axins1.bar(full_pos, vals2[0], width, color = c2, edgecolor = 'k')
axins1.bar(full_pos, vals2[1], width, bottom = vals2[0], color = c2, alpha = al, edgecolor = 'k')
if compare>2:
    axins1.bar(full_pos+width, vals3[0], width, color = c3, edgecolor = 'k')
    axins1.bar(full_pos+width, vals3[1], width, bottom = vals3[0], color = c3, alpha = al, edgecolor = 'k')
    if compare > 3:
        axins1.bar(full_pos+ 2*width, vals4[0], width, color = c4, edgecolor = 'k')
        axins1.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], color = c4, alpha = al, edgecolor = 'k')
axins1.set_xticks([4,5,6])
axins1.set_xticklabels(['DCT', 'CNT', 'CD'], fontsize=in_xlab_size)
x1 = 4 - 0.5 #CNT index
x2 = 6 + 0.75#urine index
y1 = -1
y2 = 9
axins1.set_xlim(x1,x2)
axins1.set_ylim(y1,y2)
axins1.axhline(0,color='k')

ax1.set_xticks(x+0.5*width)
ax1.set_ylabel('Na$^+$ transport ($\mu$mol/min)', fontsize=ylab_size)
ax1.set_xticklabels(segment_transport, fontsize=xlab_size)
ax1.legend(fontsize=leg_size)
ax1.text(figlab_shift, ax1.get_ylim()[1], 'A', size=fig_lab_size, weight='bold')

# #------------
# # K
# #------------
if compare == 2:
    vals1, vals2 = get_data2('K')
elif compare == 3:
    vals1, vals2, vals3 = get_data3('K')
else:
    vals1, vals2, vals3, vals4 = get_data4('K')

ax2 = fig.add_subplot(3,2,2)


# bar1
sup1 = ax2.bar(full_pos-width, vals1[0], width, align = 'center', edgecolor = 'k', color = c1, label = label1)
jux1 = ax2.bar(full_pos-width, vals1[1], width, bottom = vals1[0], align = 'center', edgecolor = 'k', color = c1, alpha = al)
# bar2
sup2 = ax2.bar(full_pos, vals2[0], width, align = 'center', edgecolor = 'k', color = c2, label = label2)
jux2 = ax2.bar(full_pos, vals2[1], width, bottom = vals2[0], align = 'center', edgecolor = 'k', color = c2, alpha = al)
if compare > 2:
    # bar 3
    sup3 = ax2.bar(full_pos+width, vals3[0], width, align = 'center', edgecolor = 'k', color = c3, label = label3)
    jux3 = ax2.bar(full_pos+width, vals3[1], width, bottom = vals3[0], align = 'center', edgecolor = 'k', color = c3, alpha = al)
    if compare > 3:
        # bar 4
        sup4 = ax2.bar(full_pos+ 2*width, vals4[0], width, align = 'center', edgecolor = 'k', color = c4, label = label4)
        jux4 = ax2.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], align = 'center', edgecolor = 'k', color = c4, alpha = al)
plt.axhline(0, color = 'k')
# inset
axins2 = inset_axes(ax2, width="30%", height = 1.1, loc='upper right')
axins2.bar(full_pos-width, vals1[0], width, color = c1, edgecolor = 'k')
axins2.bar(full_pos-width, vals1[1], width, bottom = vals1[0], color = c1, alpha = al, edgecolor = 'k')
axins2.bar(full_pos, vals2[0], width, color = c2, edgecolor = 'k')
axins2.bar(full_pos, vals2[1], width, bottom = vals2[0], color = c2, alpha = al, edgecolor = 'k')
if compare > 2:
    axins2.bar(full_pos+width, vals3[0], width, color = c3, edgecolor = 'k')
    axins2.bar(full_pos+width, vals3[1], width, bottom = vals3[0], color = c3, alpha = al, edgecolor = 'k')
    if compare > 3:
        axins2.bar(full_pos+ 2*width, vals4[0], width, color = c4, edgecolor = 'k')
        axins2.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], color = c4, alpha = al, edgecolor = 'k')
axins2.set_xticks([4,5,6])
axins2.set_xticklabels(['DCT', 'CNT', 'CD'], fontsize=in_xlab_size)
x1 = 4 - 0.5 #CNT index
x2 = 6 + 0.75#urine index
y1 = -5.0
y2 = 2.0
axins2.set_xlim(x1,x2)
axins2.set_ylim(y1,y2)
axins2.axhline(0,color='k')

ax2.set_xticks(x+0.5*width)
ax2.set_ylabel('K$^+$ transport ($\mu$mol/min)', fontsize=ylab_size)
ax2.set_xticklabels(segment_transport, fontsize=xlab_size)
#ax2.legend(fontsize=leg_size)
ax2.text(figlab_shift, ax2.get_ylim()[1], 'B', size=fig_lab_size, weight='bold')

# #------------
# # Cl
# #------------
if compare == 2:
    vals1, vals2 = get_data2('Cl')
elif compare == 3:
    vals1, vals2, vals3 = get_data3('Cl')
else:
    vals1, vals2, vals3, vals4 = get_data4('Cl')

ax3 = fig.add_subplot(3,2,3)


# bar1
sup1 = ax3.bar(full_pos-width, vals1[0], width, align = 'center', edgecolor = 'k', color = c1, label = label1)
jux1 = ax3.bar(full_pos-width, vals1[1], width, bottom = vals1[0], align = 'center', edgecolor = 'k', color = c1, alpha = al)
# bar2
sup2 = ax3.bar(full_pos, vals2[0], width, align = 'center', edgecolor = 'k', color = c2, label = label2)
jux2 = ax3.bar(full_pos, vals2[1], width, bottom = vals2[0], align = 'center', edgecolor = 'k', color = c2, alpha = al)
if compare > 2:
    # bar 3
    sup3 = ax3.bar(full_pos+width, vals3[0], width, align = 'center', edgecolor = 'k', color = c3, label = label3)
    jux3 = ax3.bar(full_pos+width, vals3[1], width, bottom = vals3[0], align = 'center', edgecolor = 'k', color = c3, alpha = al)
    if compare > 3:
        # bar 4
        sup4 = ax3.bar(full_pos+ 2*width, vals4[0], width, align = 'center', edgecolor = 'k', color = c4, label = label4)
        jux4 = ax3.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], align = 'center', edgecolor = 'k', color = c4, alpha = al)
plt.axhline(0, color = 'k')
# inset
axins3 = inset_axes(ax3, width="35%", height = 1.3, loc=5)
axins3.bar(full_pos-width, vals1[0], width, color = c1, edgecolor = 'k')
axins3.bar(full_pos-width, vals1[1], width, bottom = vals1[0], color = c1, alpha = al, edgecolor = 'k')
axins3.bar(full_pos, vals2[0], width, color = c2, edgecolor = 'k')
axins3.bar(full_pos, vals2[1], width, bottom = vals2[0], color = c2, alpha = al, edgecolor = 'k')
if compare > 2:
    axins3.bar(full_pos+width, vals3[0], width, color = c3, edgecolor = 'k')
    axins3.bar(full_pos+width, vals3[1], width, bottom = vals3[0], color = c3, alpha = al, edgecolor = 'k')
    if compare > 3:
        axins3.bar(full_pos+ 2*width, vals4[0], width, color = c4, edgecolor = 'k')
        axins3.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], color = c4, alpha = al, edgecolor = 'k')
axins3.set_xticks([4,5,6])
axins3.set_xticklabels(['DCT', 'CNT', 'CD'], fontsize=in_xlab_size)
x1 = 4 - 0.5 #CNT index
x2 = 6 + 0.75#urine index
y1 = -1
y2 = 6
axins3.set_xlim(x1,x2)
axins3.set_ylim(y1,y2)
axins3.axhline(0,color='k')

ax3.set_xticks(x+0.5*width)
ax3.set_ylabel('Cl$^-$ transport ($\mu$mol/min)', fontsize=ylab_size)
ax3.set_xticklabels(segment_transport, fontsize=xlab_size)
#ax3.legend(fontsize=leg_size)
ax3.text(figlab_shift, ax3.get_ylim()[1], 'C', size=fig_lab_size, weight='bold')

# #------------
# #NH4
# #------------
if compare == 2:
    vals1, vals2 = get_data2('NH4')
elif compare == 3:
    vals1, vals2, vals3 = get_data3('NH4')
else:
    vals1, vals2, vals3, vals4 = get_data4('NH4')

ax4 = fig.add_subplot(3,2,4)


# bar1
sup1 = ax4.bar(full_pos-width, vals1[0], width, align = 'center', edgecolor = 'k', color = c1, label = label1)
jux1 = ax4.bar(full_pos-width, vals1[1], width, bottom = vals1[0], align = 'center', edgecolor = 'k', color = c1, alpha = al)
# bar2
sup2 = ax4.bar(full_pos, vals2[0], width, align = 'center', edgecolor = 'k', color = c2, label = label2)
jux2 = ax4.bar(full_pos, vals2[1], width, bottom = vals2[0], align = 'center', edgecolor = 'k', color = c2, alpha = al)
if compare>2:
    # bar 3
    sup3 = ax4.bar(full_pos+width, vals3[0], width, align = 'center', edgecolor = 'k', color = c3, label = label3)
    jux3 = ax4.bar(full_pos+width, vals3[1], width, bottom = vals3[0], align = 'center', edgecolor = 'k', color = c3, alpha = al)
    if compare>3:
        # bar 4
        sup4 = ax4.bar(full_pos+ 2*width, vals4[0], width, align = 'center', edgecolor = 'k', color = c4, label = label4)
        jux4 = ax4.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], align = 'center', edgecolor = 'k', color = c4, alpha = al)
plt.axhline(0, color = 'k')

ax4.set_xticks(x+0.5*width)
ax4.set_ylabel('NH$_4^+$ transport ($\mu$mol/min)', fontsize=ylab_size)
ax4.set_xticklabels(segment_transport, fontsize=xlab_size)
#ax4.legend(fontsize=leg_size)
ax4.text(figlab_shift, ax4.get_ylim()[1], 'D', size=fig_lab_size, weight='bold')

# #------------
# # HCO3
# #------------
if compare == 2:
    vals1, vals2 = get_data2('HCO3')
elif compare == 3:
    vals1, vals2, vals3 = get_data3('HCO3')
else:
    vals1, vals2, vals3, vals4 = get_data4('HCO3')

ax5 = fig.add_subplot(3,2,5)


# bar1
sup1 = ax5.bar(full_pos-width, vals1[0], width, align = 'center', edgecolor = 'k', color = c1, label = label1)
jux1 = ax5.bar(full_pos-width, vals1[1], width, bottom = vals1[0], align = 'center', edgecolor = 'k', color = c1, alpha = al)
# bar2
sup2 = ax5.bar(full_pos, vals2[0], width, align = 'center', edgecolor = 'k', color = c2, label = label2)
jux2 = ax5.bar(full_pos, vals2[1], width, bottom = vals2[0], align = 'center', edgecolor = 'k', color = c2, alpha = al)
if compare>2:
    # bar 3
    sup3 = ax5.bar(full_pos+width, vals3[0], width, align = 'center', edgecolor = 'k', color = c3, label = label3)
    jux3 = ax5.bar(full_pos+width, vals3[1], width, bottom = vals3[0], align = 'center', edgecolor = 'k', color = c3, alpha = al)
    if compare>3:
        # bar 4
        sup4 = ax5.bar(full_pos+ 2*width, vals4[0], width, align = 'center', edgecolor = 'k', color = c4, label = label4)
        jux4 = ax5.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], align = 'center', edgecolor = 'k', color = c4, alpha = al)
plt.axhline(0, color = 'k')
# inset
axins5 = inset_axes(ax5, width="35%", height = 1.3, loc=5)
axins5.bar(full_pos-width, vals1[0], width, color = c1, edgecolor = 'k')
axins5.bar(full_pos-width, vals1[1], width, bottom = vals1[0], color = c1, alpha = al, edgecolor = 'k')
axins5.bar(full_pos, vals2[0], width, color = c2, edgecolor = 'k')
axins5.bar(full_pos, vals2[1], width, bottom = vals2[0], color = c2, alpha = al, edgecolor = 'k')
if compare>2:
    axins5.bar(full_pos+width, vals3[0], width, color = c3, edgecolor = 'k')
    axins5.bar(full_pos+width, vals3[1], width, bottom = vals3[0], color = c3, alpha = al, edgecolor = 'k')
    if compare>3:
        axins5.bar(full_pos+ 2*width, vals4[0], width, color = c4, edgecolor = 'k')
        axins5.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], color = c4, alpha = al, edgecolor = 'k')
axins5.set_xticks([4,5,6])
axins5.set_xticklabels(['DCT', 'CNT', 'CD'], fontsize=in_xlab_size)
x1 = 4 - 0.5 #CNT index
x2 = 6 + 0.75#urine index
y1 = -0.2
y2 = 0.75
axins5.set_xlim(x1,x2)
axins5.set_ylim(y1,y2)
axins5.axhline(0,color='k')

ax5.set_xticks(x)
ax5.set_ylabel('HCO$_3^-$ transport ($\mu$mol/min)', fontsize=ylab_size)
ax5.set_xticklabels(segment_transport, fontsize=xlab_size)
#ax5.legend(fontsize=leg_size)
ax5.text(figlab_shift, ax5.get_ylim()[1], 'E', size=fig_lab_size, weight='bold')

# #------------
# # Volume
# #------------
if compare==2:
    vals1, vals2 = get_data2('Volume')
elif compare == 3:
    vals1, vals2, vals3 = get_data3('Volume')
else:
    vals1, vals2, vals3, vals4 = get_data4('Volume')

ax6 = fig.add_subplot(3,2,6)


# bar1
sup1 = ax6.bar(full_pos-width, vals1[0], width, align = 'center', edgecolor = 'k', color = c1, label = label1)
jux1 = ax6.bar(full_pos-width, vals1[1], width, bottom = vals1[0], align = 'center', edgecolor = 'k', color = c1, alpha = al)
# bar2
sup2 = ax6.bar(full_pos, vals2[0], width, align = 'center', edgecolor = 'k', color = c2, label = label2)
jux2 = ax6.bar(full_pos, vals2[1], width, bottom = vals2[0], align = 'center', edgecolor = 'k', color = c2, alpha = al)
if compare>2:
    # bar 3
    sup3 = ax6.bar(full_pos+width, vals3[0], width, align = 'center', edgecolor = 'k', color = c3, label = label3)
    jux3 = ax6.bar(full_pos+width, vals3[1], width, bottom = vals3[0], align = 'center', edgecolor = 'k', color = c3, alpha = al)
    if compare>3:
        # bar 4
        sup4 = ax6.bar(full_pos+ 2*width, vals4[0], width, align = 'center', edgecolor = 'k', color = c4, label = label4)
        jux4 = ax6.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], align = 'center', edgecolor = 'k', color = c4, alpha = al)
plt.axhline(0, color = 'k')
# inset
axins6 = inset_axes(ax6, width="35%", height = 1.3, loc=5)
axins6.bar(full_pos-width, vals1[0], width, color = c1, edgecolor = 'k')
axins6.bar(full_pos-width, vals1[1], width, bottom = vals1[0], color = c1, alpha = al, edgecolor = 'k')
axins6.bar(full_pos, vals2[0], width, color = c2, edgecolor = 'k')
axins6.bar(full_pos, vals2[1], width, bottom = vals2[0], color = c2, alpha = al, edgecolor = 'k')
if compare>2:
    axins6.bar(full_pos+width, vals3[0], width, color = c3, edgecolor = 'k')
    axins6.bar(full_pos+width, vals3[1], width, bottom = vals3[0], color = c3, alpha = al, edgecolor = 'k')
    if compare>3:
        axins6.bar(full_pos+ 2*width, vals4[0], width, color = c4, edgecolor = 'k')
        axins6.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], color = c4, alpha = al, edgecolor = 'k')
axins6.set_xticks([4,5,6])
axins6.set_xticklabels(['DCT', 'CNT', 'CD'], fontsize=in_xlab_size)
x1 = 4 - 0.5 #CNT index
x2 = 6 + 0.75#urine index
y1 = -0.01
y2 = 0.15
axins6.set_xlim(x1,x2)
axins6.set_ylim(y1,y2)
axins6.axhline(0,color='k')

ax6.set_xticks(x)
ax6.set_ylabel('Volume transport (ml/min)', fontsize=ylab_size)
ax6.set_xticklabels(segment_transport, fontsize=xlab_size)
#ax6.legend(fontsize=leg_size)
ax6.text(figlab_shift, ax6.get_ylim()[1], 'F', size=fig_lab_size, weight='bold')






