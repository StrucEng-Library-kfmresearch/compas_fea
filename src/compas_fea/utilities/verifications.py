# Author(s): Sophia V. Kuhn (ETHZ), Marius  Weber (ETHZ, HSLU T&A)

#maybe used in the future otherwise TODO remove
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from compas.datastructures import Mesh
from compas.geometry import distance_point_point
from compas.topology import dijkstra_path
from compas.utilities import geometric_key

from time import time




try:
    import numpy as np
except ImportError:
    pass


try:
    # from scipy.interpolate import griddata
    from scipy.sparse import csr_matrix
except ImportError:
    pass


# actually used
import statistics
import math


def verification(structure, step, field='all', D_max=None, tau_cd=None):
    """Extracts analysis results from structure.results and calculated verifications and saves these in the results.
    Parameters
    ----------
    structure : obj
        The Structure object to read from.
    step : str
        name of the load step to take the analysis results from (e.g. 'step_1').
    fields : str
        verification field requests (e.g. 'all', 'shear').
    D_max: float
        maximum aggregate size [mm] (e.g. 32)
    tau_cd: float
        Design value of shear stress limit [N/mm2] (e.g. 1.4)

    Returns
    -------
    None.
    """

    # ==============================================================================
    # shear verification
    # ==============================================================================

    if field == 'shear' or field == 'all':

        if structure.results[step]=='ERROR': #if analysis has failed (no convergence)
            print('Due to the analysis error no shear verification can not be done. ')

        else:
            if D_max != None and tau_cd != None: #required imput atm

                #initialize #TODO when analysis fails then ?Error is printed in results then no shear verification can be done! solve!
                structure.results[step]['element']['v_rd']={}
                structure.results[step]['element']['v0']={}
                structure.results[step]['element']['eta_v']={}
                structure.results[step]['element']['eps_06d_loc_decisive']={}
                structure.results[step]['element']['dv']={}
                structure.results[step]['element']['centroid']={}

                # extract data
                kmax = structure.element_count() # Anzahl Elemente, Startwert bei 1 nicht bei 0!
                #data = structure.results[step]['element_info']
                
        
                
                # Extract data for SHEAR VERIFICATION
                data_eps_x_06d_bot = structure.results[step]['GP']['eps_x_06d_bot']   
                data_eps_y_06d_bot = structure.results[step]['GP']['eps_y_06d_bot']   
                data_eps_xy_06d_bot = structure.results[step]['GP']['eps_xy_06d_bot']   
                
                data_eps_x_06d_top = structure.results[step]['GP']['eps_x_06d_top']   
                data_eps_y_06d_top = structure.results[step]['GP']['eps_y_06d_top']   
                data_eps_xy_06d_top = structure.results[step]['GP']['eps_xy_06d_top']   

                # psi_1=structure.results[step]['element_info']['psi_1'].values() 
                # psi_2=structure.results[step]['element_info']['psi_2'].values() 
                # psi_3=structure.results[step]['element_info']['psi_3'].values() 
                # psi_4=structure.results[step]['element_info']['psi_4'].values() 

                dm_1=structure.results[step]['element_info']['dm1'].values() 
                dm_2=structure.results[step]['element_info']['dm2'].values() 
                dm_3=structure.results[step]['element_info']['dm3'].values() 
                dm_4=structure.results[step]['element_info']['dm4'].values()   

                oo=structure.results[step]['element_info']['oo'].values() 
                uu=structure.results[step]['element_info']['uu'].values() 
                h_shell=structure.results[step]['element_info']['h_shell'].values() 
                LNr_d06_bot=structure.results[step]['element_info']['LNr_d06_bot'].values() 
                LNr_d06_top=structure.results[step]['element_info']['LNr_d06_top'].values()      
                # fcc=structure.results[step]['element_info']['fcc'].values() 

                # Start 
                for k in range(kmax): 
                    ele_type = structure.results[step]['element']['ele_type'][k].values()

                    #initialize
                    v_rd=None
                    v0=None

                    if ele_type[0] == 1.0:
                        #sf4-->vx
                        sf4 = structure.results[step]['element']['sf4'][k].values()
                        vx = statistics.mean(list(sf4)) # Bildet den durchschnitt aller Integrationspunkte des Elements
                        #sf5-->vy
                        sf5 = structure.results[step]['element']['sf5'][k].values()
                        vy = statistics.mean(list(sf5)) # Bildet den durchschnitt aller Integrationspunkte des Elements
                        centroid= structure.element_centroid(element=k)   
                        #vx+vy-->v0
                        v0 = (vx**2+vy**2)**0.5
                        
                        angle_v0=math.atan(vy/vx) # Vorzeichen zur Berechnung der Darstellung gegenuber dem Usermat umkehren

                        # step 5: Shear strength and verification (optional)
                        anlge_shear_ver=angle_v0 # Independon from a "main reinforcement"

                        # Calculate mean values fro eps_06d bot and top (original values on GP and not at element midpoint!)    
                        k_eps=k*4   
                        eps_x_06d_bot_mean=(data_eps_x_06d_bot[k_eps]+data_eps_x_06d_bot[k_eps+1]+data_eps_x_06d_bot[k_eps+2]+data_eps_x_06d_bot[k_eps+3])/4
                        eps_y_06d_bot_mean=(data_eps_y_06d_bot[k_eps]+data_eps_y_06d_bot[k_eps+1]+data_eps_y_06d_bot[k_eps+2]+data_eps_y_06d_bot[k_eps+3])/4
                        eps_xy_06d_bot_mean=(data_eps_xy_06d_bot[k_eps]+data_eps_xy_06d_bot[k_eps+1]+data_eps_xy_06d_bot[k_eps+2]+data_eps_xy_06d_bot[k_eps+3])/4
                        eps_06d_bot_mean=eps_x_06d_bot_mean*math.cos(anlge_shear_ver)**2+eps_y_06d_bot_mean*math.sin(anlge_shear_ver)**2+eps_xy_06d_bot_mean*math.cos(anlge_shear_ver)*math.sin(anlge_shear_ver)

                        eps_x_06d_top_mean=(data_eps_x_06d_top[k_eps]+data_eps_x_06d_top[k_eps+1]+data_eps_x_06d_top[k_eps+2]+data_eps_x_06d_top[k_eps+3])/4
                        eps_y_06d_top_mean=(data_eps_y_06d_top[k_eps]+data_eps_y_06d_top[k_eps+1]+data_eps_y_06d_top[k_eps+2]+data_eps_y_06d_top[k_eps+3])/4
                        eps_xy_06d_top_mean=(data_eps_xy_06d_top[k_eps]+data_eps_xy_06d_top[k_eps+1]+data_eps_xy_06d_top[k_eps+2]+data_eps_xy_06d_top[k_eps+3])/4
                        eps_06d_top_mean=eps_x_06d_top_mean*math.cos(anlge_shear_ver)**2+eps_y_06d_top_mean*math.sin(anlge_shear_ver)**2+eps_xy_06d_top_mean*math.cos(anlge_shear_ver)*math.sin(anlge_shear_ver)
                        
                        # Decision which values (top or bot) is decisive and calculation of the corresponding d_v (midpoint element)
                        if eps_06d_top_mean > 0 and eps_06d_bot_mean < 0:
                            #case: only top is positive ->eps from top
                            eps_06d_mean=eps_06d_top_mean 
                            eps_06d_loc_decisive='top'
                            d_s1_top = h_shell[k]-uu[k]-dm_1[k]/2 
                            d_s2_top = h_shell[k]-uu[k]-dm_1[k]-dm_2[k]/2
                            d_v=(d_s1_top+d_s2_top)/2 # from top
                            Nr_layer=LNr_d06_top[k]
                        elif eps_06d_top_mean > 0 and eps_06d_bot_mean > 0:
                            #case: both bot and top are positive ->take maximum eps
                            if eps_06d_top_mean >= eps_06d_bot_mean:
                                eps_06d_mean=eps_06d_top_mean  
                                eps_06d_loc_decisive='top' 
                                d_s1_top = h_shell[k]-uu[k]-dm_1[k]/2 
                                d_s2_top = h_shell[k]-uu[k]-dm_1[k]-dm_2[k]/2
                                d_v=(d_s1_top+d_s2_top)/2 # from top    
                                Nr_layer=LNr_d06_top[k]        
                            elif eps_06d_top_mean < eps_06d_bot_mean:
                                eps_06d_mean=eps_06d_bot_mean 
                                eps_06d_loc_decisive='bot'
                                d_s4_bot = h_shell[k]-oo[k]-dm_4[k]/2 
                                d_s3_bot = h_shell[k]-oo[k]-dm_4[k]-dm_3[k]/2
                                d_v=(d_s4_bot+d_s3_bot)/2 # from bot 
                                Nr_layer=LNr_d06_bot[k]                   
                        elif eps_06d_top_mean <= 0 and eps_06d_bot_mean <= 0:
                            #case: both bot and top are negative ->take no eps
                            eps_06d_mean=None #not defined as positiv
                            eps_06d_loc_decisive=None
                            d_v=None
                        elif eps_06d_top_mean < 0 and eps_06d_bot_mean > 0:
                            #case: only bot is positive ->take eps from bot
                            eps_06d_mean=eps_06d_bot_mean
                            eps_06d_loc_decisive='bot'
                            d_s4_bot = h_shell[k]-oo[k]-dm_4[k]/2 
                            d_s3_bot = h_shell[k]-oo[k]-dm_4[k]-dm_3[k]/2
                            d_v=(d_s4_bot+d_s3_bot)/2 # from bot  
                            Nr_layer=LNr_d06_bot[k]

                        # Calculation of shear resistance accorindg to SIA 269              
                        if eps_06d_mean != None:    
                            k_g=48/(16+D_max)
                            k_d=1/(1+2.5*eps_06d_mean*d_v*k_g)                
                            v_rd=k_d*tau_cd*d_v                    
                        else:
                            v_rd=None

                        #save verification results of element k to results dict
                        structure.results[step]['element']['v_rd'][k]=v_rd
                        structure.results[step]['element']['v0'][k]=v0
                        if v_rd != None and v0 !=None:
                            structure.results[step]['element']['eta_v'][k]=v_rd/v0
                        else:
                            structure.results[step]['element']['eta_v'][k]=None
                        structure.results[step]['element']['eps_06d_loc_decisive'][k]=eps_06d_loc_decisive
                        structure.results[step]['element']['dv'][k]=d_v
                        structure.results[step]['element']['centroid'][k]= structure.element_centroid(element=k) 
                    


            else: 
                raise ValueError('D_max or tau_cd is None. Please define D_max and tau_cd. This is necessary for the shear verficiation')
        

    