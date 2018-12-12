#############################################################################################
#Documentation: In neuroscience, fiber photometry combined with bulk calcium imaging is     #
#a method utilized to record the activity of a defined neuronal population in a target brain# 
#area from a free-moving animal performing a given task. In this method, a virus containing #
# a calcium indicator (i.e. AAV-DIO-GCaMP) is injected in the brain region of interest and  #
# an optical fiber is implanted nearby the brain region of interest that now expresses the  #
#calcium indicator. Through the optical fiber, an excitation light is delivered and the     # 
#overall fluorescence induced by calcium-activity during the light excitation is collected. #
#Fiber photometry sums up the overall fluorescence of neurons expressing the genetically    #
#encoded calcium indicator, GCaMP, and can be used as a proxy for neuronal activation.      #
#The typical setup for freely behaving animals consists of the excitation light source      #
# (470 nm - to excite GCaMP), the beamsplitter/combiner that separates excitation and       #
#fluorescence light, the fiber-optic rotary joint, the optical cannula and connecting       #
#fiber-optic patch cords. Our set-up also contains an autofluorescence light (405 nm) which #
#is used as a control light to account for photobleaching and movement artifacts.           #
#For more info, please see Beas et al 2018.                                                 #
#																							#
#The purpose of this code is to analyze and calculate ΔF/F0 for photometry recordings from  #
#a mouse expressing GCaMP in the paraventricular nucleus of the thalamus (PVT)              #
#as in (Beas et al 2018). The activity of neurons of the PVT were recorded while the mouse  #
#ate food pellets.																			#
#																							#
#After imaging from mice, three CSV files were exported from the photometry set up:			#
# (These files were uploaded on supplemental material as a sample files for testing 		#
# the code)																					#
# 1- autofluorescence signal ("A1_food2_AF.csv").											#
# 2- GCaMP signal ("A1_food2_GC.csv").														#
#		Note:(our photometry rig collects 8 samples per second for both the GCaMP and 		#
# 		the autofluorescence) 																#
# 3- TTL pulses showing the times when the mouse ate a food pellet ("A1_food2_pulse.csv")	#
#																							#
#ΔF/F0 analysis: As in Beas et al 2008, this code will analyze data by first applying  		# 
#a least-squares linear fit to the 405-nm signal (auto) to align it to the 470-nm 			#
#signal (GCaMP).  																			#
#The resulting fitted 405-nm signal was then used to normalize the 473-nm signal 			#
#as follows:																				#
# ΔF/F0 = (473-nm signal − fitted 405-nm signal)/fitted 405-nm signal. 						#
#																							#
#############################################################################################



# Import packages needed
import pandas as pd
import numpy as np
from pylab import *
import matplotlib.lines as mlines #matplotlib
import matplotlib.pyplot as plt #matplotlib
plt.style.use('ggplot') #emulate ggplot from R
from sklearn import preprocessing

# Locating the times in the imaging ("GCaMP and Auto") when the mouse ate a pellet to calculate ΔF/F0 
def locate_pulses(pt, master):
    results = []
    for element in pt:
        subset = master[master['time'] > element]
        record = subset[subset['time'] < element + 0.2]
        results.append(record['time'].values)
    return results

# Create a method for ΔF/F0 and label it as "DF_F"
def DF_F_analysis(indx, lower, upper, m_df):

    # create pulses_onset
    pulse_onset = indx  
    if pulse_onset in pulse_onset_rows:
        print('pulse_onset =', pulse_onset)
    else:
        raise ValueError('pulse_onset not found in pulse_onset_rows')
    # create begin_input
    begin_input = pulse_onset - lower  
    print('begin_input =', begin_input)
    # create last_row
    last_row = len(m_df) - pulse_onset - 1
    # create end_input
    end_input = pulse_onset + upper  
    print('end_input =', end_input)
    # create file name for future files
    file_name = 'A1_pellet_' + str(indx)  # IMPORTANT - to change file name, format as 'file name'
    print(file_name)

    # create master_input pandas DataFrame
    master_input = m_df[begin_input:end_input + 1]

    # create master_trendlines pandas DataFrame
    master_trendlines = master_input.loc[begin_input:pulse_onset - 1]
    # reset master_trendlines row index
    master_trendlines = master_trendlines.reset_index(drop=True)
    # create x_master_trendlines (ranging from 1 to # rows in master_trendlines) pandas DataFrame
    x_range_master_trendlines = master_trendlines.axes[0] - (master_trendlines.axes[0][0] - 1)
    x_master_trendlines = pd.DataFrame({'x': x_range_master_trendlines})
    # add x_master_trendlines to master_trendlines
    master_trendlines = pd.concat([x_master_trendlines, master_trendlines], axis=1, join='inner')

    # determine gcamp linear trendline equation
    (a, b) = polyfit(master_trendlines.x, master_trendlines.auto_d0, 1)
    auto_linear_trendline_equation = 'y = ' + str(round(a, 5)) + 'x + ' + str(round(b, 5))

    # create auto scatter plot with auto linear trendline
    plt.scatter(master_trendlines.x, master_trendlines.auto_d0, color='blue', s=10)
    auto_trendline_values = polyval([a, b], master_trendlines.x)
    plt.plot(master_trendlines.x, auto_trendline_values, linewidth=3, color='red')
    plt.title(auto_linear_trendline_equation, fontsize=10, y=0.9)
    plt.xlabel('x')
    plt.ylabel('autofluorescence')
    # save auto scatter plot as PDF to working directory
    plt.show()
    plt.savefig('auto_plot_' + file_name + '.pdf')

    # determine gcamp linear trendline equation
    (c, d) = polyfit(master_trendlines.x, master_trendlines.gc_d0, 1)
    gcamp_linear_trendline_equation = 'y = ' + str(round(c, 5)) + 'x + ' + str(round(d, 5))

    # create gcamp scatter plot with gcamp linear trendline
    plt.scatter(master_trendlines.x, master_trendlines.gc_d0, color='blue', s=10)
    gcamp_trendline_values = polyval([c, d], master_trendlines.x)
    plt.plot(master_trendlines.x, gcamp_trendline_values, linewidth=3, color='red')
    plt.title(gcamp_linear_trendline_equation, fontsize=10, y=0.9)
    plt.xlabel('x')
    plt.ylabel('gcamp')
    # save gcamp scatter plot as PDF to working directory
    plt.show()
    plt.savefig('gcamp_plot_' + file_name + '.pdf')

    # create master_calculations pandas DataFrame
    master_calculations = master_input
    # reset master_calculations row index
    master_calculations = master_calculations.reset_index(drop=True)
    # create x_master_calculations (ranging from 1 to # rows in master_calculations) pandas DataFrame
    x_range_master_calculations = master_calculations.axes[0] - (master_calculations.axes[0][0] - 1)
    x_master_calculations = pd.DataFrame({'x': x_range_master_calculations})
    # add x_master_calculations to master_calculations
    master_calculations = pd.concat([x_master_calculations, master_calculations], axis=1, join='inner')

    # create auto_trendline_y pandas DataFrame
    auto_y_list = []
    for x in x_range_master_calculations:
        y = a * x + b
        auto_y_list.append(y)
    auto_y_array = np.array(auto_y_list)
    auto_trendline_y = pd.DataFrame({'auto_trendline_y': auto_y_array})
    # add auto_trendline_y to master_calculations
    master_calculations = pd.concat([master_calculations, auto_trendline_y], axis=1, join='inner')
    # create gcamp_trendline_y pandas DataFrame
    gcamp_y_list = []
    for x in x_range_master_calculations:
        y = c * x + d
        gcamp_y_list.append(y)
    gcamp_y_array = np.array(gcamp_y_list)
    gcamp_trendline_y = pd.DataFrame({'gcamp_trendline_y': gcamp_y_array})
    # add gcamp_trendline_y to master_calculations
    master_calculations = pd.concat([master_calculations, gcamp_trendline_y], axis=1, join='inner')

    # subtract auto_trendline_y from auto to create auto_fit column in master_calculations
    master_calculations['auto_fit'] = master_calculations['auto_d0'] - master_calculations['auto_trendline_y']
    # subtract gcamp_trendline_y from gcamp to create gcamp_fit column in master_calculations
    master_calculations['gcamp_fit'] = master_calculations['gc_d0'] - master_calculations['gcamp_trendline_y']

    # add gcamp_trendline y-intercept to auto_fit to create auto_fit column in master_calculations
    master_calculations['auto_final'] = d + master_calculations['auto_fit']
    # add gcamp_trendline y-intercept to gcamp_fit to create gcamp_fit column in master_calculations
    master_calculations['gcamp_final'] = d + master_calculations['gcamp_fit']

    # calculate delta f/f (dff) and create dff column in master_calculations
    master_calculations['dff'] = ((master_calculations['gcamp_final'] - master_calculations['auto_final']) /
                                  master_calculations['auto_final']) * 100

    # write out master_calculations as a csv file to working directory
    master_calculations.to_csv('master_calculations_' + file_name + '.csv')

    # determine the rows in which pulse occurs
    pulse_rows_dff = master_calculations.loc[master_calculations.pulses_d0 > 1].index[:].tolist()
    # determine the rows in which pulse onset occurs
    pulse_onset_rows_dff = [pulse_rows_dff[0]]
    for i in range(2, len(pulse_rows_dff)):
        if pulse_rows_dff[i] > pulse_rows_dff[i - 1] + 1:
            pulse_onset_rows_dff.append(pulse_rows_dff[i])
    # determine the times in which pulse onset occurs
    pulse_onset_time_dff = list(master_calculations.time.loc[pulse_onset_rows_dff])
    # create dff line plot
    fig, ax = plt.subplots()
    ax.plot(master_calculations.time, master_calculations.dff, color='blue')
    ax.set_xlabel('time (sec)')
    ax.set_ylabel('delta f/f')
    x = 0
    while x < len(pulse_onset_time_dff):
        ax.annotate(' ', xy=(pulse_onset_time_dff[x], min(master_calculations.dff)),
                    arrowprops=dict(facecolor='black', shrink=0.05))
        x = x + 1
    arrow = mlines.Line2D([], [], color='black', marker='^', markersize=12, label='pulse onset')
    ax.legend(handles=[arrow])
    # save dff line plot as PDF to working directory
    fig.savefig('dff_plot_' + file_name + '.pdf')

    print('Success')

if __name__ == '__main__':

    auto = pd.read_csv('./A1_food2_AF.csv') #Import autofluorescence signal file
    gc = pd.read_csv('./A1_food2_GC.csv') #Import GCaMP signal file
    pulses = pd.read_csv('./A1_food2_pulse.csv') #Import times (as pulses) when mouse ate a food pellet

    # Slice data frames into only columns of interest.
    time = auto['TIME'].values # Time is needed to be able to determine when the mouse ate the pellet.
    auto_d0 = auto['D0'].values # Raw autofluorescence signal
    gc_d0 = gc['D0'].values # Raw GCaMP signal
    master_df = pd.DataFrame(columns=('time', 'auto_d0', 'gc_d0'), index=np.arange(0,len(time)))

    # Creating master data frame containing only the important columns needed for DF_F analysis from both auto and gcamp
    master_df['time'] = time
    master_df['auto_d0'] = auto_d0
    master_df['gc_d0'] = gc_d0

    # Insert new column (pulses_d0) into master dataframe. Initialize with zeros.
    pulse_d0 = np.zeros(len(time))

    # Insert pulses_d0 into master dataframe
    master_df['pulses_d0'] = pulse_d0

    pulses_time = pulses['TIME'].values

    #round time to be able to localizate it better in our pulse time (since photometry rig takes 8 samples per second)
    pulses_time = [round(t, 1) for t in pulses_time]


    #checking point for loacalizing the rows where pulses are located
    sub_set = master_df[master_df['time'] > 136]

    p_results = locate_pulses(pulses_time, master_df)

    #Getting rid of duplicate values found
    buffer = []
    for el in p_results:
        buffer.append(el[0])

    p_results = buffer

    # Mark down in pulses column with a 5 the time when the pulse occurred (replace the 0 with a 5).
    for el in p_results:
        indx = master_df[master_df['time'] == el].index[0]
        master_df.loc[indx]['pulses_d0'] = 5

    master_df.to_csv('./master_with_pulses.csv')

    plot_df = pd.DataFrame(master_df)
    plot_df = plot_df.set_index('time')
    plot_df = plot_df.drop('auto_d0', axis=1)


    #Scale gcamp signal for better visualization
    plot_df['gc_d0'] = preprocessing.scale(plot_df['gc_d0'].values)

    fig = plt.figure(figsize=(10, 5))
    ax = fig.gca()
    for indx, col in enumerate(plot_df.columns):
        plt.plot(plot_df[plot_df.index > 2][col], label=col)

    plt.title('GCamp Raw')
    plt.ylabel("Scaled GC_D0", fontsize=14)
    plt.xticks(rotation=45)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('gcamp_raw_A1.pdf')

    # Checking that the 0 in pulses_d0 was replaced with a 5 in particular rows where the pulse happened
    subset = master_df[master_df['pulses_d0'] == 5]
    subset.to_csv('./A1_pulse_onset_rows.csv')

    pulse_onset_rows = [i for i in subset.index[:]]


    # Showing where (what rows) the rest of the pulses are located
    subset = master_df[master_df['pulses_d0'] == 5]
    print(subset)

	#Use method to calculate DF_F for each time the mouse ate a food pellet
    DF_F_analysis(indx=1081, lower=40, upper=80, m_df=master_df) #Important: Change indx value for each "pulse". 
    #Change lower to establish a baseline (pre-pulse). Change upper to establish the length of the analysis (Post-pulse). 
    #Reminder: Our photometry rig takes 8 samples per second. Thus, for a 5 second baseline we need a lower of "40"  
    #and for a 10 sec post pulse analysis, we need an upper of "80".
    
    
    
    
    