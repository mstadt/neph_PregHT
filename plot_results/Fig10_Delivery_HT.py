
import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#--------------------------------------------------------------------
# begin user input
#--------------------------------------------------------------------
compare = 4 # 2, 3, or 4

solute_list = ['Na', 'K', 'Cl', 'HCO3', 'NH4', 'Volume']

segs_early = ['PT', 'SDL', 'mTAL', 'cTAL', 'DCT', 'CNT']
segs_cd = ['CCD', 'IMCD']
seg_labels = ['PT', 'DL', 'mTAL', 'cTAL', 'DCT', 'CNT', 'CCD', 'urine']
# 8 segments


direct1 = 'NP-baseline2022-04-22'
sex1 = 'female'

direct2 = 'NP-HTN'
sex2 = 'female'

direct3 = 'midpregnant_rat'
sex3 = 'midpregnant'

direct4 = 'MP-HTN'
sex4 = 'midpregnant'

label1 = 'Normotensive virgin'
label2 = 'Hypertensive virgin'
label3 = 'Normotensive MP'
label4 = 'Hypertensive MP'

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


#-----------------------------------------------------------------------
#end user input
#----------------------------------------------------------------------
#===========================================================================
# functions used
#===========================================================================

def get_delivery(direct, sex, solute, segment, supOrjux):
    os.chdir(direct)
    if solute == 'Volume':
        fname = sex+'_'+humOrrat+'_'+segment+'_water_volume_in_Lumen'+supOrjux+'.txt'
        file = open(fname, 'r')
        delivery = float(file.readline())
        file.close()
    else:
        fname = sex+'_'+humOrrat+'_'+segment+'_flow_of_'+solute+'_in_Lumen'+supOrjux+'.txt'
        file = open(fname, 'r')
        delivery = float(file.readline())
        file.close()
    os.chdir('..')
    # convert to per kidney in micromoles
    delivery_cf = delivery * cf
    return delivery_cf

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

def get_cd_data(direct, sex, solute, segments):
    # this is for the CD
    direct_deliv = []
    for seg in segments:
        if seg[-2:].lower() != 'cd':
            raise Exception('only for the collecting duct')
        direct_deliv.append(float(get_delivery(direct, sex, solute, seg, '')))
    direct_deliv_cf = direct_deliv * cf
    return direct_deliv_cf

def get_data(direct, sex, solute, segments):
    direct_deliv_sup = []
    direct_deliv_jux1 = []
    direct_deliv_jux2 = []
    direct_deliv_jux3 = []
    direct_deliv_jux4 = []
    direct_deliv_jux5 = []
    
    for seg in segments:
        direct_deliv_sup.append(float(get_delivery(direct, sex, solute, seg, '_sup')))
        direct_deliv_jux1.append(float(get_delivery(direct, sex, solute, seg, '_jux1')))
        direct_deliv_jux2.append(float(get_delivery(direct, sex, solute, seg, '_jux2')))
        direct_deliv_jux3.append(float(get_delivery(direct, sex, solute, seg, '_jux3')))
        direct_deliv_jux4.append(float(get_delivery(direct, sex, solute, seg, '_jux4')))
        direct_deliv_jux5.append(float(get_delivery(direct, sex, solute, seg, '_jux5')))
    
        if seg.lower() == 'cnt':
            direct_deliv_sup.append(float(get_out_deliv(direct, sex, solute, 'cnt', '_sup')))
            direct_deliv_jux1.append(float(get_out_deliv(direct,sex, solute, 'cnt', '_jux1')))
            direct_deliv_jux2.append(float(get_out_deliv(direct, sex, solute, 'cnt', '_jux2')))
            direct_deliv_jux3.append(float(get_out_deliv(direct, sex, solute, 'cnt', '_jux3')))
            direct_deliv_jux4.append(float(get_out_deliv(direct, sex, solute, 'cnt', '_jux4')))
            direct_deliv_jux5.append(float(get_out_deliv(direct, sex, solute, 'cnt', '_jux5')))
        
    
    temp= np.array(direct_deliv_sup)*neph_weight[0] + np.array(direct_deliv_jux1)*neph_weight[1]\
        + np.array(direct_deliv_jux2)*neph_weight[2] + np.array(direct_deliv_jux3)*neph_weight[3] \
            + np.array(direct_deliv_jux4)*neph_weight[4] + np.array(direct_deliv_jux4)*neph_weight[5]
    direct_deliv_number = temp
    
    # row 0 is superficial vals
    # row 1 is jux1, row 2 jux2,...,row5 jux5
    direct_vals = np.matrix([direct_deliv_sup, direct_deliv_jux1, direct_deliv_jux2, 
                            direct_deliv_jux3, direct_deliv_jux4, direct_deliv_jux5])
    
    sup_temp = direct_vals[0] * neph_weight[0]
    sup_temp1 = np.array(sup_temp)
    sup_vals_weight = sup_temp1[0]

    jux_temp = direct_vals[1]*neph_weight[1] + direct_vals[2]*neph_weight[2]+\
        direct_vals[3]*neph_weight[3] + direct_vals[4]*neph_weight[4] + direct_vals[5]*neph_weight[5]
    jux_temp1 = np.array(jux_temp)
    jux_vals_weight = jux_temp1[0]
    
    return direct_deliv_number, direct_vals, sup_vals_weight, jux_vals_weight

def get_vals2(solute):
    direct_deliv_num1, direct_vals1, sup_vals1, jux_vals1 = get_data(direct1, sex1, solute, segs_early)
    direct_deliv_num2, direct_vals2, sup_vals2, jux_vals2 = get_data(direct2, sex2, solute, segs_early)
    
    urine1 = get_out_deliv(direct1, sex1, solute, 'IMCD','')
    urine2 = get_out_deliv(direct2, sex2, solute, 'IMCD', '')
    
    
    vals1 = np.zeros(2*8).reshape(2,8)
    vals2 = np.zeros(2*8).reshape(2,8)
    for i in range(7):
        vals1[0][i] = sup_vals1[i]
        vals1[1][i] = jux_vals1[i]
        
        vals2[0][i] = sup_vals2[i]
        vals2[1][i] = jux_vals2[i]
        
    vals1[0][7] = urine1
    vals2[0][7] = urine2

    
    return vals1, vals2

def get_vals3(solute):
    direct_deliv_num1, direct_vals1, sup_vals1, jux_vals1 = get_data(direct1, sex1, solute, segs_early)
    direct_deliv_num2, direct_vals2, sup_vals2, jux_vals2 = get_data(direct2, sex2, solute, segs_early)
    direct_deliv_num3, direct_vals3, sup_vals3, jux_vals3 = get_data(direct3, sex3, solute, segs_early)
    
    urine1 = get_out_deliv(direct1, sex1, solute, 'IMCD','')
    urine2 = get_out_deliv(direct2, sex2, solute, 'IMCD', '')
    urine3 = get_out_deliv(direct3, sex3, solute, 'IMCD', '')
    
    
    vals1 = np.zeros(2*8).reshape(2,8)
    vals2 = np.zeros(2*8).reshape(2,8)
    vals3 = np.zeros(2*8).reshape(2,8)

    for i in range(7):
        vals1[0][i] = sup_vals1[i]
        vals1[1][i] = jux_vals1[i]
        
        vals2[0][i] = sup_vals2[i]
        vals2[1][i] = jux_vals2[i]
        
        vals3[0][i] = sup_vals3[i]
        vals3[1][i] = jux_vals3[i]

    vals1[0][7] = urine1
    vals2[0][7] = urine2
    vals3[0][7] = urine3

    
    return vals1, vals2, vals3


def get_vals4(solute):
    direct_deliv_num1, direct_vals1, sup_vals1, jux_vals1 = get_data(direct1, sex1, solute, segs_early)
    direct_deliv_num2, direct_vals2, sup_vals2, jux_vals2 = get_data(direct2, sex2, solute, segs_early)
    direct_deliv_num3, direct_vals3, sup_vals3, jux_vals3 = get_data(direct3, sex3, solute, segs_early)
    direct_deliv_num4, direct_vals4, sup_vals4, jux_vals4 = get_data(direct4, sex4, solute, segs_early)
    
    urine1 = get_out_deliv(direct1, sex1, solute, 'IMCD','')
    urine2 = get_out_deliv(direct2, sex2, solute, 'IMCD', '')
    urine3 = get_out_deliv(direct3, sex3, solute, 'IMCD', '')
    urine4 = get_out_deliv(direct4, sex4, solute, 'IMCD', '')
    
    
    vals1 = np.zeros(2*8).reshape(2,8)
    vals2 = np.zeros(2*8).reshape(2,8)
    vals3 = np.zeros(2*8).reshape(2,8)
    vals4 = np.zeros(2*8).reshape(2,8)
    for i in range(7):
        vals1[0][i] = sup_vals1[i]
        vals1[1][i] = jux_vals1[i]
        
        vals2[0][i] = sup_vals2[i]
        vals2[1][i] = jux_vals2[i]
        
        vals3[0][i] = sup_vals3[i]
        vals3[1][i] = jux_vals3[i]
        
        vals4[0][i] = sup_vals4[i]
        vals4[1][i] = jux_vals4[i]
    vals1[0][7] = urine1
    vals2[0][7] = urine2
    vals3[0][7] = urine3
    vals4[0][7] = urine4
    
    return vals1, vals2, vals3, vals4
                    
#===============================================================================
#   make the figures
#===============================================================================
# colors
c1 = 'c'
c2 = 'darkcyan'
c3 = 'mediumvioletred'
c4 = 'palevioletred'
al = 0.4

title_size = 25
ylab_size = 20
leg_size=18
fig_lab_size = 25
xlab_size=20
in_xlab_size=16

# figure specs
figlab_shift = -1.5
width = 0.2

figwidth = 15
figheight = 20
fig = plt.figure(figsize=(figwidth,figheight))
plt.rcParams.update({'font.size':20})


# positiions
sup_pos = np.arange(len(seg_labels[:7]))
full_pos = np.arange(len(seg_labels))
later_pos = np.arange(len(seg_labels[:7]), len(seg_labels))

#------------
# Na
#------------
if compare == 2:
    vals1, vals2 = get_vals2('Na')
elif compare == 3:
    vals1, vals2, vals3 = get_vals3('Na')
else:
    vals1, vals2, vals3, vals4 = get_vals4('Na')




ax1 = fig.add_subplot(3,1,1)


# bar1
sup1 = ax1.bar(full_pos-width, vals1[0], width, align = 'center', edgecolor = 'k', color = c1, label = label1)
jux1 = ax1.bar(full_pos-width, vals1[1], width, bottom = vals1[0], align = 'center', edgecolor = 'k', color = 'white')
# bar2
sup2 = ax1.bar(full_pos, vals2[0], width, align = 'center', edgecolor = 'k', color = c2, label = label2)
jux2 = ax1.bar(full_pos, vals2[1], width, bottom = vals2[0], align = 'center', edgecolor = 'k', color = 'white')
if compare>2:
    # bar 3
    sup3 = ax1.bar(full_pos+width, vals3[0], width, align = 'center', edgecolor = 'k', color = c3, label = label3)
    jux3 = ax1.bar(full_pos+width, vals3[1], width, bottom = vals3[0], align = 'center', edgecolor = 'k', color = 'white')
    if compare>3:
        # bar 4
        sup4 = ax1.bar(full_pos+ 2*width, vals4[0], width, align = 'center', edgecolor = 'k', color = c4, label = label4)
        jux4 = ax1.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], align = 'center', edgecolor = 'k', color = 'white')
# inset
axins1 = inset_axes(ax1, width="30%", height = 1.5, loc=5)
axins1.bar(full_pos-width, vals1[0], width, color = c1, edgecolor = 'k')
axins1.bar(full_pos-width, vals1[1], width, bottom = vals1[0], edgecolor = 'k', color = 'white')
axins1.bar(full_pos, vals2[0], width, color = c2, edgecolor = 'k')
axins1.bar(full_pos, vals2[1], width, bottom = vals2[0], edgecolor = 'k', color = 'white')
if compare>2:
    axins1.bar(full_pos+width, vals3[0], width, color = c3, edgecolor = 'k')
    axins1.bar(full_pos+width, vals3[1], width, bottom = vals3[0], edgecolor = 'k', color = 'white')
    if compare>3:
        axins1.bar(full_pos+ 2*width, vals4[0], width, color = c4, edgecolor = 'k')
        axins1.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0],  edgecolor = 'k', color = 'white')
axins1.set_xticks(np.arange(4,8) + 0.5*width)
axins1.set_xticklabels(['DCT', 'CNT', 'CCD', 'urine'], fontsize=in_xlab_size)
x1 = 4 - 2*width #CNT index
x2 = 7 + 3*width #urine index
y1 = 0
y2 = 12
axins1.set_xlim(x1,x2)
axins1.set_ylim(y1,y2)

ax1.set_xticks(full_pos+0.5*width)
ax1.set_ylabel('Na$^+$ delivery ($\mu$mol/min)', fontsize=ylab_size)
ax1.set_xticklabels(seg_labels, fontsize=xlab_size)
ax1.legend(fontsize=leg_size, loc = (0.25,0.6))
ax1.text(figlab_shift, ax1.get_ylim()[1], 'A', size=fig_lab_size, weight='bold')

#------------
# K
#------------
if compare == 2:
    vals1, vals2 = get_vals2('K')
elif compare == 3:
    vals1, vals2, vals3 = get_vals3('K')
else:
    vals1, vals2, vals3, vals4 = get_vals4('K')

ax2 = fig.add_subplot(3,1,2)


# bar1
sup1 = ax2.bar(full_pos-width, vals1[0], width, align = 'center', edgecolor = 'k', color = c1, label = label1)
jux1 = ax2.bar(full_pos-width, vals1[1], width, bottom = vals1[0], align = 'center', edgecolor = 'k', color = 'white')
# bar2
sup2 = ax2.bar(full_pos, vals2[0], width, align = 'center', edgecolor = 'k', color = c2, label = label2)
jux2 = ax2.bar(full_pos, vals2[1], width, bottom = vals2[0], align = 'center', edgecolor = 'k', color = 'white')
if compare>2:
    # bar 3
    sup3 = ax2.bar(full_pos+width, vals3[0], width, align = 'center', edgecolor = 'k', color = c3, label = label3)
    jux3 = ax2.bar(full_pos+width, vals3[1], width, bottom = vals3[0], align = 'center', edgecolor = 'k', color = 'white')
    if compare>3:
        # bar 4
        sup4 = ax2.bar(full_pos+ 2*width, vals4[0], width, align = 'center', edgecolor = 'k', color = c4, label = label4)
        jux4 = ax2.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], align = 'center', edgecolor = 'k', color = 'white')

ax2.set_xticks(full_pos+0.5*width)
ax2.set_ylabel('K$^+$ delivery ($\mu$mol/min)', fontsize=ylab_size)
ax2.set_xticklabels(seg_labels, fontsize=xlab_size)
ax2.legend(fontsize=leg_size, loc = (0.25,0.6))
ax2.text(figlab_shift, ax2.get_ylim()[1], 'B', size=fig_lab_size, weight='bold')



# #------------
# # Volume
# #------------
if compare == 2:
    vals1, vals2 = get_vals2('Volume')
elif compare == 3:
    vals1, vals2, vals3 = get_vals3('Volume')
elif compare == 4:
    vals1, vals2, vals3, vals4 = get_vals4('Volume')

ax6 = fig.add_subplot(3,1,3)


# bar1
sup1 = ax6.bar(full_pos-width, vals1[0], width, align = 'center', edgecolor = 'k', color = c1, label = label1)
jux1 = ax6.bar(full_pos-width, vals1[1], width, bottom = vals1[0], align = 'center', edgecolor = 'k', color = 'white')
# bar2
sup2 = ax6.bar(full_pos, vals2[0], width, align = 'center', edgecolor = 'k', color = c2, label = label2)
jux2 = ax6.bar(full_pos, vals2[1], width, bottom = vals2[0], align = 'center', edgecolor = 'k', color = 'white')
if compare > 2:
    # bar 3
    sup3 = ax6.bar(full_pos+width, vals3[0], width, align = 'center', edgecolor = 'k', color = c3, label = label3)
    jux3 = ax6.bar(full_pos+width, vals3[1], width, bottom = vals3[0], align = 'center', edgecolor = 'k', color = 'white')
    if compare > 3:
        # bar 4
        sup4 = ax6.bar(full_pos+ 2*width, vals4[0], width, align = 'center', edgecolor = 'k', color = c4, label = label4)
        jux4 = ax6.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], align = 'center', edgecolor = 'k', color = 'white')

# inset
axins6 = inset_axes(ax6, width="30%", height = 1.5, loc=5)
axins6.bar(full_pos-width, vals1[0], width, color = c1, edgecolor = 'k')
axins6.bar(full_pos-width, vals1[1], width, bottom = vals1[0], color = 'white', edgecolor = 'k')
axins6.bar(full_pos, vals2[0], width, color = c2, edgecolor = 'k')
axins6.bar(full_pos, vals2[1], width, bottom = vals2[0], color = 'white', edgecolor = 'k')
if compare > 2:
    axins6.bar(full_pos+width, vals3[0], width, color = c3, edgecolor = 'k')
    axins6.bar(full_pos+width, vals3[1], width, bottom = vals3[0], color = 'white', edgecolor = 'k')
    if compare > 3:
        axins6.bar(full_pos+ 2*width, vals4[0], width, color = c4, edgecolor = 'k')
        axins6.bar(full_pos+ 2*width, vals4[1], width, bottom = vals4[0], color = 'white', edgecolor = 'k')
axins6.set_xticks(np.arange(4,8) + 0.5*width)
axins6.set_xticklabels(['DCT', 'CNT', 'CCD', 'urine'], fontsize=in_xlab_size)
x1 = 4 - 2*width #CNT index
x2 = 7 + 3*width #urine index
y1 = 0
y2 = 0.275
axins6.set_xlim(x1,x2)
axins6.set_ylim(y1,y2)


ax6.set_xticks(full_pos+0.5*width)
ax6.set_ylabel('Volume delivery (ml/min)', fontsize=ylab_size)
ax6.set_xticklabels(seg_labels, fontsize=xlab_size)
ax6.legend(fontsize=leg_size, loc = (0.25,0.6))
ax6.text(figlab_shift, ax6.get_ylim()[1], 'C', size=fig_lab_size, weight='bold')