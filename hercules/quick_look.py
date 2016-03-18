#name: quick_look.py
#author: Will Barnes
#date created: 8 March 2016

#description: Class for quickly and easily creating plots of HYDRAD data

#import needed modules
import os,sys
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import seaborn.apionly as sns

class QuickViewer(object):
    """
    Class for quickly and easily examining HYDRAD results
    
    Attributes
    ----------
    hydrad_results_root: string
        root directory where relevant HYDRAD results are stored
    
    Examples
    --------
    
    Notes
    --------
    * make_profile_movie() function heavily based on this blog post 
     (https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/) from Jake VanderPlas
    
    """
    
    def __init__(self,hydrad_results_root,dpi=1000, fontsize=18., figsize=(8,8), fontsize_eps=0.75, fformat='pdf', **kwargs):
        self.hydrad_results_root = hydrad_results_root
        #set plot styling options
        self.dpi,self.fontsize,self.figsize,self.fontsize_eps,self.fformat = dpi,fontsize,figsize,fontsize_eps,fformat
        #configure logger
        self.logger = logging.getLogger(type(self).__name__)
        
    
    def load_results(self):
        """Read in HYDRAD results"""
        
        #count number of relevant files
        self.num_phy_files = len([pf for pf in os.listdir(self.hydrad_results_root) if '.phy' in pf])
        if self.num_phy_files == 0:
            raise ValueError('No *.phy files found in %s'%self.hydrad_results_root)
        #loop over files
        self.results = []
        self.time = []
        for i in range(self.num_phy_files):
            #load parameters
            tmp = np.loadtxt(os.path.join(self.hydrad_results_root,'profile%d.phy'%i))
            self.results.append({'s':tmp[:,0], 'Te':tmp[:,-4], 'Ti':tmp[:,-3], 'ne':tmp[:,3], 'ni':tmp[:,4], 'pe':tmp[:,5],'pi':tmp[:,6],'v':tmp[:,1]})
            #load time information
            with open(os.path.join(self.hydrad_results_root,'profile%d.amr'%i),'r') as f:
                self.time.append(float(f.readline()))
                
        #convert to numpy array
        self.time = np.array(self.time)
            
            
    def plot_profile(self,indices=None,line_cm='coolwarm',print_fig_filename=None,**kwargs):
        """Plot snapshots in time of loop profiles"""
        
        if not hasattr(self,'results'):
            raise AttributeError('No results found. Make sure to run self.load_results() first.')
        
        colors = sns.color_palette(line_cm,self.num_phy_files)
        
        if indices is None:
            indices = range(self.num_phy_files)
            
        #setup figure
        fig,axes = plt.subplots(2,2,figsize=self.figsize,sharex=True)
         
        #plotting   
        for i in indices:
            r = self.results[i]
            #temperature
            axes[0,0].plot(r['s'],r['Te']/1e+6,color=colors[i],linestyle='-',alpha=0.5)
            axes[0,0].plot(r['s'],r['Ti']/1e+6,color=colors[i],linestyle='--',alpha=0.5)
            #density
            axes[0,1].plot(r['s'],r['ne'],color=colors[i],linestyle='-',alpha=0.5)
            axes[0,1].plot(r['s'],r['ni'],color=colors[i],linestyle='--',alpha=0.5)
            #velocity
            axes[1,0].plot(r['s'],r['v'],color=colors[i],linestyle='-',alpha=0.5)
            #pressure
            axes[1,1].plot(r['s'],r['pe'],color=colors[i],linestyle='-',alpha=0.5)
            axes[1,1].plot(r['s'],r['pi'],color=colors[i],linestyle='--',alpha=0.5)
            
        #plot styling
        axes[0,0].set_ylabel(r'$T$ $(\mathrm{MK})$',fontsize=self.fontsize)
        axes[0,1].set_ylabel(r'$n$  $(\mathrm{cm}^{-3})$',fontsize=self.fontsize)
        axes[1,0].set_ylabel(r'$v$ $(\mathrm{cm/s})$',fontsize=self.fontsize)
        axes[1,1].set_ylabel(r'$p$ $(\mathrm{dyne}\,\mathrm{cm}^{-2})$',fontsize=self.fontsize)
        axes[1,0].set_xlabel(r'$s$ $(\mathrm{cm})$',fontsize=self.fontsize)
        axes[1,1].set_xlabel(r'$s$ $(\mathrm{cm})$',fontsize=self.fontsize)
        axes[0,1].set_yscale('log')
        axes[1,1].set_yscale('log')
        axes[0,0].set_xlim([r['s'][0],r['s'][-1]])
        
        plt.tight_layout()
        
        #save or show the figure
        if print_fig_filename is not None:
            plt.savefig(print_fig_filename+'.'+self.fformat,format=self.fformat,dpi=self.dpi)
        else:
            plt.show()
            
            
    def make_profile_movie(self,line_cm='coolwarm',movie_filename=None,interval=20,fps=30,**kwargs):
        """Make movie showing loop profile at each time t"""
        
        if not hasattr(self,'results'):
            raise AttributeError('No results found. Make sure to run self.load_results() first.')
        
        #colors
        colors = sns.color_palette(line_cm,self.num_phy_files)
        
        #setup figure
        fig,axes = plt.subplots(2,2,figsize=self.figsize,sharex=True)
        #find and set static axes limits
        min_t,max_t=1.e+300,0.
        min_n,max_n=1.e+300,0.
        min_v,max_v=1.e+300,0.
        min_p,max_p=1.e+300,0.
        for r in self.results:
            #temperature
            if r['Te'].min() < min_t or r['Ti'].min() < min_t:
                min_t=min(r['Te'].min(),r['Ti'].min())
            if r['Te'].max() > max_t or r['Ti'].max() > max_t:
                max_t=max(r['Te'].max(),r['Ti'].max())
            #density
            if r['ne'].min() < min_n or r['ni'].min() < min_n:
                min_n=min(r['ne'].min(),r['ni'].min())
            if r['ne'].max() > max_n or r['ni'].max() > max_n:
                max_n=max(r['ne'].max(),r['ni'].max())
            #velocity
            if r['v'].min() < min_v:
                min_v=r['v'].min()
            if r['v'].max() > max_v:
                max_v=r['v'].max() 
            #pressure
            if r['pe'].min() < min_p or r['pi'].min() < min_p:
                min_p=min(r['pe'].min(),r['pi'].min())
            if r['pe'].max() > max_p or r['pi'].max() > max_p:
                max_p=max(r['pe'].max(),r['pi'].max())
            
        #y axes
        axes[0,0].set_ylim([min_t,max_t])
        axes[0,1].set_ylim([min_n,max_n])
        axes[1,0].set_ylim([min_v,max_v])
        axes[1,1].set_ylim([min_p,max_p])        
        #x axes    
        axes[0,0].set_xlim([self.results[0]['s'][0],self.results[0]['s'][-1]])
        axes[0,1].set_xlim([self.results[0]['s'][0],self.results[0]['s'][-1]])
        axes[1,0].set_xlim([self.results[0]['s'][0],self.results[0]['s'][-1]])
        axes[1,1].set_xlim([self.results[0]['s'][0],self.results[0]['s'][-1]])
        
        #labels
        axes[0,0].set_ylabel(r'$T$ $(\mathrm{MK})$',fontsize=self.fontsize)
        axes[0,1].set_ylabel(r'$n$  $(\mathrm{cm}^{-3})$',fontsize=self.fontsize)
        axes[1,0].set_ylabel(r'$v$ $(\mathrm{cm/s})$',fontsize=self.fontsize)
        axes[1,1].set_ylabel(r'$p$ $(\mathrm{dyne}\,\mathrm{cm}^{-2})$',fontsize=self.fontsize)
        axes[1,0].set_xlabel(r'$s$ $(\mathrm{cm})$',fontsize=self.fontsize)
        axes[1,1].set_xlabel(r'$s$ $(\mathrm{cm})$',fontsize=self.fontsize)
        axes[0,1].set_yscale('log')
        axes[1,1].set_yscale('log')
        
        plt.tight_layout()
        
        #setup plotters for each axes
        lineTe,=axes[0,0].plot([],[],linestyle='-',linewidth=2)
        lineTi,=axes[0,0].plot([],[],linestyle='--',linewidth=2)
        linene,=axes[0,1].plot([],[],linestyle='-',linewidth=2)
        lineni,=axes[0,1].plot([],[],linestyle='--',linewidth=2)
        linev,=axes[1,0].plot([],[],linestyle='-',linewidth=2)
        linepe,=axes[1,1].plot([],[],linestyle='-',linewidth=2)
        linepi,=axes[1,1].plot([],[],linestyle='--',linewidth=2)
        
        #set title with current time
        time_title=fig.suptitle('',fontsize=self.fontsize)
        
        #initialization function
        def init():
            """initialize animation"""
            time_title.set_text('')
            #temperature
            lineTe.set_data([])
            lineTe.set_color([])
            lineTi.set_data([])
            lineTi.set_color([])
            #density
            linene.set_data([])
            linene.set_color([])
            lineni.set_data([])
            lineni.set_color([])
            #velocity
            linev.set_data([],[])
            linev.set_color([])
            #pressure
            linepe.set_data([])
            linepe.set_color([])
            linepi.set_data([])
            linepi.set_color([])
        
        #animation function
        def animate(i):
            r=self.results[i]
            time_title.set_text(r'$t=%.2f$'%(self.time[i]))
            #temperature
            lineTe.set_data(r['s'],r['Te']/1e+6)
            lineTe.set_color(colors[i])
            lineTi.set_data(r['s'],r['Ti']/1e+6)
            lineTi.set_color(colors[i])
            #density
            linene.set_data(r['s'],r['ne'])
            linene.set_color(colors[i])
            lineni.set_data(r['s'],r['ni'])
            lineni.set_color(colors[i])
            #velocity
            linev.set_data(r['s'],r['v'])
            linev.set_color(colors[i])
            #pressure
            linepe.set_data(r['s'],r['pe'])
            linepe.set_color(colors[i])
            linepi.set_data(r['s'],r['pi'])
            linepi.set_color(colors[i])
            return lineTe,lineTi,linene,lineni,linev,linepe,linepi,time_title
            
        #make the animation
        anim=animation.FuncAnimation(fig,animate,init_func=init,frames=self.num_phy_files,blit=True,interval=interval,)
        
        #save or show
        if movie_filename is not None:
            anim.save(movie_filename,fps=fps,extra_args=['-vcodec','libx264'])
        else:
            plt.show()
            
            
    def make_timeseries(self,lower_percent=0.25,upper_percent=0.25,**kwargs):
        """Average properties over loop length"""
        
        #read in loop length
        with open(os.path.join(self.hydrad_results_root,'profile0.amr'),'r') as f:
            _=f.readline()
            _=f.readline()
            self.L=float(f.readline())
            
        #get bounds    
        lower_bound = self.L/2. - lower_percent*self.L/2.
        upper_bound = self.L/2. + upper_percent*self.L/2.
        
        self.timeseries = {'Te':[],'Ti':[],'ne':[],'ni':[],'v':[],'pe':[],'pi':[]}
        #iterate over results
        for r in self.results:
            #find relevant indices
            ind_avg = np.where((r['s']>=lower_bound) & (r['s']<=upper_bound))[0]
            if len(ind_avg) <= 0:
                self.logger.error("Cannot find average for bounds=(%f,%f)"%(lower_bound,upper_bound))
                continue
            s_avg = r['s'][ind_avg]
            for key in r:
                if key != 's':
                    self.timeseries[key].append(np.average(r[key][ind_avg], weights=np.gradient(s_avg)))
                    
        #convert to numpy arrays
        for key in self.timeseries:
            self.timeseries[key] = np.array(self.timeseries[key])

            
    def plot_timeseries(self,print_fig_filename=None,**kwargs):
        """Plot average parameters as a function of time"""
        
        if not hasattr(self,'timeseries'):
            raise AttributeError('No timeseries found. Make sure to run self.make_timeseries() first.')
        
        #create figure
        fig,axes = plt.subplots(3,1,figsize=self.figsize,sharex=True)
        #plot data
        axes[0].plot(self.time ,self.timeseries['Te']/1e+6, color=sns.color_palette('deep')[0], linestyle='-',linewidth=2)
        axes[0].plot(self.time, self.timeseries['Ti']/1e+6, color=sns.color_palette('deep')[2], linestyle='-',linewidth=2)
        axes[1].plot(self.time ,self.timeseries['ne']/1e+8, color=sns.color_palette('deep')[0], linestyle='-',linewidth=2)
        axes[1].plot(self.time, self.timeseries['ni']/1e+8, color=sns.color_palette('deep')[2], linestyle='-',linewidth=2)
        axes[2].plot(self.time ,self.timeseries['pe'], color=sns.color_palette('deep')[0], linestyle='-',linewidth=2)
        axes[2].plot(self.time, self.timeseries['pi'], color=sns.color_palette('deep')[2], linestyle='-',linewidth=2)
        
        #plot styling
        axes[0].set_ylabel(r'$T$ $(\mathrm{MK})$',fontsize=self.fontsize)
        axes[1].set_ylabel(r'$n$  $(10^8\,\mathrm{cm}^{-3})$',fontsize=self.fontsize)
        axes[2].set_ylabel(r'$p$ $(\mathrm{dyne}\,\mathrm{cm}^{-2})$',fontsize=self.fontsize)
        axes[2].set_xlabel(r'$t$ $(\mathrm{s})$',fontsize=self.fontsize)
        axes[0].set_xlim([self.time[0],self.time[-1]])
        
        plt.tight_layout()
        
        #save or show the figure
        if print_fig_filename is not None:
            plt.savefig(print_fig_filename+'.'+self.fformat,format=self.fformat,dpi=self.dpi)
        else:
            plt.show()
