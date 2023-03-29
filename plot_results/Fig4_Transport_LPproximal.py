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

compare = 3

direct1 = 'latepregnant_rat'
sex1 = 'latepregnant'

direct2 = 'LP-PTlen'
sex2 = 'latepregnant'

direct3 = 'LP-NHE3'
sex3 = 'latepregnant'





label1 = 'Baseline LP'
label2 = 'LP PT length at virgin'
label3 = 'LP NHE3 at virgin'




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

def get_out_deliv(direct, sex, solute, segment, supOrjux):
    # flow at the end of given segment
    os.chdir(direct)
    if solute == 'Volume':
        fname = sex+'_'+humOrrat+'_'+segment+'_water_volume_in_Lumen'+supOrjux+'.txt'
        file = open(fname, 'r')
        delivery = float(np.loadtxt(fname, delimiter = '\n', unpack = True)[-1])
        file.close()
    else:
        fname = sex+'_'+humOrrat+'_'+segment+'_flow_of_'+solute+'_in_Lumen'+supOrjux+'.txt'
        file = open(fname, 'r')
        delivery = float(np.loadtxt(fname, delimiter = '\n', unpack = True)[-1])
        file.close()
    os.chdir('..')
    # convert
    delivery_cf = cf * delivery
    return delivery_cf




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
     
    
    
    urine1 = get_out_deliv(direct1, sex1, solute, 'IMCD','')
    urine2 = get_out_deliv(direct2, sex2, solute, 'IMCD', '')
    urine3 = get_out_deliv(direct3, sex3, solute, 'IMCD', '')

    
    DTvals1 = np.zeros(2*8).reshape(2,8)
    DTvals1[0][0:7] = sup1_final
    DTvals1[0][-1] = urine1
    DTvals1[1][0:7] = jux1_final
    
    DTvals2 = np.zeros(2*8).reshape(2,8)
    DTvals2[0][0:7] = sup2_final
    DTvals2[0][-1] = urine2
    DTvals2[1][0:7] = jux2_final
    
    DTvals3 = np.zeros(2*8).reshape(2,8)
    DTvals3[0][0:7] = sup3_final
    DTvals3[0][-1] = urine3
    DTvals3[1][0:7] = jux3_final
    

    

    
    return DTvals1, DTvals2, DTvals3
    
                    
#===============================================================================
#   make the figures
#===============================================================================
# colors
c1 = 'green'
c2 = 'springgreen'
c3 = 'yellow'
c4 = 'darkgreen'

juxc = 'white'
al = 0.4

title_size = 25
ylab_size = 20
leg_size=18
fig_lab_size = 25
xlab_size=20
in_xlab_size=16

# figure specs
figlab_shift = -1.75
width = 0.2


figheight = 15
figwidth  = 10
fig = plt.figure(figsize = (figwidth,figheight))
plt.rcParams.update({'font.size':20})


x = np.arange(8)
full_pos = np.arange(8)

xlabels = ['PT', 'DL', 'LAL', 'TAL', 'DCT', 'CNT', 'CD', 'urine']

#------------
# Na
#------------

vals1, vals2, vals3 = get_data3('Na')


ax1 = fig.add_subplot(2,1,1)

# bar1
sup1 = ax1.bar(full_pos-width, vals1[0], width, align = 'center', edgecolor = 'k', color = c1, label = label1)
jux1 = ax1.bar(full_pos-width, vals1[1], width, bottom = vals1[0], align = 'center', edgecolor = 'k', color = 'white')
# bar2
sup2 = ax1.bar(full_pos, vals2[0], width, align = 'center', edgecolor = 'k', color = c2, label = label2)
jux2 = ax1.bar(full_pos, vals2[1], width, bottom = vals2[0], align = 'center', edgecolor = 'k', color = 'white')

# bar 3
sup3 = ax1.bar(full_pos+width, vals3[0], width, align = 'center', edgecolor = 'k', color = c3, label = label3)
temp = vals3[0]
temp[1] = 0 
jux3 = ax1.bar(full_pos+width, vals3[1], width, bottom = temp, align = 'center', edgecolor = 'k', color = 'white')



plt.axhline(0, color = 'k')
# inset
axins1 = inset_axes(ax1, width="35%", height = 1.0, loc = 5)
axins1.bar(full_pos-width, vals1[0], width, color = c1, edgecolor = 'k')
axins1.bar(full_pos-width, vals1[1], width, bottom = vals1[0], color ='white', edgecolor = 'k')
axins1.bar(full_pos, vals2[0], width, color = c2, edgecolor = 'k')
axins1.bar(full_pos, vals2[1], width, bottom = vals2[0], color = 'white', edgecolor = 'k')
if compare>2:
    axins1.bar(full_pos+width, vals3[0], width, color = c3, edgecolor = 'k')
    axins1.bar(full_pos+width, vals3[1], width, bottom = vals3[0], color = 'white', edgecolor = 'k')
axins1.set_xticks([4,5,6,7])
axins1.set_xticklabels(['DCT', 'CNT', 'CD', 'urine'], fontsize=in_xlab_size)
x1 = 4 - 0.5 #CNT index
x2 = 7 + 0.5#urine index
y1 = -1
y2 = 6
axins1.set_xlim(x1,x2)
axins1.set_ylim(y1,y2)
axins1.axhline(0,color='k')

ax1.set_xticks(x)
ax1.set_ylabel('Na$^+$ transport ($\mu$mol/min)', fontsize=ylab_size)
ax1.set_xticklabels(xlabels, fontsize=xlab_size)
ax1.legend(fontsize=leg_size)
ax1.set_ylim(-5.0,105)
ax1.text(figlab_shift, ax1.get_ylim()[1], 'A', size=fig_lab_size, weight='bold')

# #------------
# # K
# #------------
vals1, vals2, vals3= get_data3('K')

ax2 = fig.add_subplot(2,1,2)


# bar1
sup1 = ax2.bar(full_pos-width, vals1[0], width, align = 'center', edgecolor = 'k', color = c1, label = label1)
jux1 = ax2.bar(full_pos-width, vals1[1], width, bottom = vals1[0], align = 'center', edgecolor = 'k', color = 'white')
# bar2
sup2 = ax2.bar(full_pos, vals2[0], width, align = 'center', edgecolor = 'k', color = c2, label = label2)
jux2 = ax2.bar(full_pos, vals2[1], width, bottom = vals2[0], align = 'center', edgecolor = 'k', color = 'white')

# bar 3
sup3 = ax2.bar(full_pos+width, vals3[0], width, align = 'center', edgecolor = 'k', color = c3, label = label3)
jux3 = ax2.bar(full_pos+width, vals3[1], width, bottom = vals3[0], align = 'center', edgecolor = 'k', color = 'white')



plt.axhline(0, color = 'k')
# inset
axins2 = inset_axes(ax2, width="27%", height = 1.1, loc='upper right')
axins2.bar(full_pos-width, vals1[0], width, color = c1, edgecolor = 'k')
axins2.bar(full_pos-width, vals1[1], width, bottom = vals1[0], color = 'white', edgecolor = 'k')
axins2.bar(full_pos, vals2[0], width, color = c2, edgecolor = 'k')
axins2.bar(full_pos, vals2[1], width, bottom = vals2[0], color = 'white', edgecolor = 'k')
if compare > 2:
    axins2.bar(full_pos+width, vals3[0], width, color = c3, edgecolor = 'k')
    axins2.bar(full_pos+width, vals3[1], width, bottom = vals3[0], color = 'white', edgecolor = 'k')
axins2.set_xticks([4,5,6])
axins2.set_xticklabels(['DCT', 'CNT', 'CD'], fontsize=in_xlab_size)
x1 = 4 - 0.5 #CNT index
x2 = 6 + 0.5
y1 = -0.75
y2 = 0.1
axins2.set_xlim(x1,x2)
axins2.set_ylim(y1,y2)
axins2.axhline(0,color='k')

ax2.set_xticks(x)
ax2.set_ylabel('K$^+$ transport ($\mu$mol/min)', fontsize=ylab_size)
ax2.set_xticklabels(xlabels, fontsize=xlab_size)
ax2.legend(fontsize=leg_size, loc = (0.15,0.74))
ax2.set_ylim(-1.0,3.5)
ax2.text(figlab_shift, ax2.get_ylim()[1], 'B', size=fig_lab_size, weight='bold')








