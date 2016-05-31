def interpolate_timecourse(label_model,t,tf,nt,size,condition):
   maxt=float(tf) #tf is the final time of the simulation, assuming simulation starts at t0
   nt=float(nt)   # nt is the number of time points in the simulation
   dt=(maxt)/(nt-1)     # how much time increases between time points. 1 is substracted because it includes time 0 
   n=t/dt         
   n_low=int(n) #Converts decimal to int ie 4.6 ->4
   n_up=n_low+1
   label_model.lower_size_dict[size]
   #print [n_low,n_up]
   #print [t,n_low,round(t,8),maxt]
   for local_size in label_model.lower_size_dict[size]:
       for mi in label_model.lower_size_dict[size][local_size]:
           n_variable=label_model.size_variable_dict[local_size][mi]
           if round(t,8)==0:
              label_model.condition_size_yy_dict[condition][local_size][n_variable]=label_model.condition_size_yy0_dict[condition][local_size][n_variable]
           elif round(t,8)>=maxt:
              label_model.condition_size_yy_dict[condition][local_size][n_variable]=round(label_model.size_timecoure_dict[local_size][nt-1][n_variable],8)
           else:
              begining=label_model.size_timecoure_dict[local_size][n_low][n_variable]
              end=label_model.size_timecoure_dict[local_size][n_up][n_variable]
              m=(end-begining)/dt
              local_t=dt*n_low+t
              label_model.condition_size_yy_dict[condition][local_size][n_variable]=round(begining+m*local_t,8)
           #print label_model.condition_size_yy_dict[condition][size][n_variable]"""


def interpolate_timecourse(label_model,local_t,size,condition):
    for local_size in label_model.lower_size_dict[size]:
        for mi in label_model.lower_size_dict[size][local_size]:
            n_variable=label_model.size_variable_dict[local_size][mi]
            xp = t
            fp=label_model.size_timecoure_dict[local_size][:,n_variable]
            label_model.condition_size_yy_dict[condition][local_size][n_variable]=np.interp(local_t, xp, fp)  
