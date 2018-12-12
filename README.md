# Beas_Final_Project
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
#																							                                              #
#The purpose of this code is to analyze and calculate ΔF/F0 for photometry recordings from  #
#a mouse expressing GCaMP in the paraventricular nucleus of the thalamus (PVT)              #
#as in (Beas et al 2018). The activity of neurons of the PVT were recorded while the mouse  #
#ate food pellets.																			                                    #
#																							                                              #
#After imaging from mice, three CSV files were exported from the photometry set up:			    #
# (These files were uploaded on supplemental material as a sample files for testing 		    #
# the code)																				                                        	#
# 1- autofluorescence signal ("A1_food2_AF.csv").											                      #
# 2- GCaMP signal ("A1_food2_GC.csv").													                          	#
#		Note:(our photometry rig collects 8 samples per second for both the GCaMP and 		      #
# 		the autofluorescence) 																                                #
# 3- TTL pulses showing the times when the mouse ate a food pellet ("A1_food2_pulse.csv")	  #
#																						                                              	#
#ΔF/F0 analysis: As in Beas et al 2008, this code will analyze data by first applying  		  # 
#a least-squares linear fit to the 405-nm signal (auto) to align it to the 470-nm 			    #
#signal (GCaMP).  																			                                    #
#The resulting fitted 405-nm signal was then used to normalize the 473-nm signal 			      #
#as follows:																				                                        #
# ΔF/F0 = (473-nm signal − fitted 405-nm signal)/fitted 405-nm signal. 						          #
#																						                                              	#
#############################################################################################
