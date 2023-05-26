# Author(s): Andrew Liew (github.com/andrewliew), Marius Weber (IBK, ETHZ)

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



__all__ = [
    'Materials',
]

MPa = 10**(-6)
GPa = 10**(-9)


class Materials(object):

    def __init__(self):
        pass

    def write_materials(self):
    
        self.write_section('Materials')
        self.blank_line()
        
        materials = self.structure.materials
        element_set = self.structure.sets
        properties = self.structure.element_properties
        sections = self.structure.sections
        

        # Lesen der secnum
        ele_prop=[]
        
        for key_ele_prop in sorted(properties):
            ele_prop.append(key_ele_prop)

        # Lesen der benotigten Geometrien
        t_set=[]
        nn_layer_set=[]
        for key_prop in sorted(properties):
            #self.write_subsection(key_prop)

            property = properties[key_prop]
            section = sections[property.section]
            geometry = section.geometry  
            nr_layers=section.nr_layers   

            if nr_layers is not None:                
                nn_layer_set.append(nr_layers.get('nn', None))
            else:
                nn_layer_set.append(None)

            if geometry is not None:
                t_set.append(geometry.get('t', None))
            else:
                t_set.append(None)
        
        
        # Lesen der vorhandenen element sets
        ele_set=[]
        for key_set in sorted(self.structure.sets):
                        
            element_set = self.structure.sets[key_set]
            if element_set.type != 'node': 
                #self.write_subsection(key_set)
                ele_set.append(key_set)

        # Start Schlaufe zur Zuweisung Material zu Element sets                
        count_mat=-1
        for key, material in sorted(materials.items()):
            


            # ------------------------------------------------------------------------------------------------------
            # Material and element properties einlesen
            # ------------------------------------------------------------------------------------------------------
            count_mat=count_mat+1
            mtype = material.__name__

            if mtype == 'ElasticIsotropic':  
                E = material.E['E']
                v = material.v['v']
                p = material.p   
                t=t_set[count_mat]  
                nn=nn_layer_set[count_mat] 
                G=E/(2*(1+v))
                E111=5/6*G*t    
                E221=5/6*G*t                     
                E121=0                
                delta_h=t/nn
                
                self.write_subsection(key)
            
            elif mtype == 'MPCStiff':  

                self.blank_line()
                self.write_line('!No Material properties for MPCs are needed') 
                self.blank_line()                                
            
            elif mtype == 'CMMUsermat':                  
                
                
                # Dies noch anpassen delten



                # HIER WEITRFAHREN MIT ANDEREN PARAMETER
                R_Rohr = material.R_Rohr['R_Rohr']
                rho = material.rho['rho']
                oo = material.oo['oo']
                uu = material.uu['uu']
                
                beton = material.beton['beton']
                fcc = material.fcc['fcc']                
                vc = material.vc['vc']
                ecu = material.ecu['ecu']
                k_E = material.k_E['k_E']
                theta_b0 = material.theta_b0['theta_b0']
                theta_b1 = material.theta_b1['theta_b1']
                k_riss = material.k_riss['k_riss']
                Entfestigung = material.Entfestigung['Entfestigung']
                lambdaTS = material.lambdaTS['lambdaTS']
                srmx = material.srmx['srmx']
                srmy = material.srmy['srmy']
                Begrenzung = material.Begrenzung['Begrenzung']
                KritQ = material.KritQ['KritQ']
                winkelD = material.winkelD['winkelD']
                k_vr = material.k_vr['k_vr']
                fswy = material.fswy['fswy']


                stahl1 = material.stahl1['stahl1']
                zm1 = material.zm1['zm1']
                fsy1 = material.fsy1['fsy1']
                fsu1 = material.fsu1['fsu1']
                esu1 = material.esu1['esu1']
                esv1 = material.esv1['esv1']
                Es1 = material.Es1['Es1']
                ka1 = material.ka1['ka1']
                kb1 = material.kb1['kb1']
                kc1 = material.kc1['kc1']
                as1 = material.as1['as1']
                dm1 = material.dm1['dm1']
                psi1 = material.psi1['psi1']

                stahl2 = material.stahl2['stahl2']
                zm2 = material.zm2['zm2']
                fsy2 = material.fsy2['fsy2']
                fsu2 = material.fsu2['fsu2']
                esu2 = material.esu2['esu2']
                esv2 = material.esv2['esv2']
                Es2 = material.Es2['Es2']
                ka2 = material.ka2['ka2']
                kb2 = material.kb2['kb2']
                kc2 = material.kc2['kc2']
                as2 = material.as2['as2']
                dm2 = material.dm2['dm2']
                psi2 = material.psi2['psi2']

                stahl3 = material.stahl3['stahl3']
                zm3 = material.zm3['zm3']
                fsy3 = material.fsy3['fsy3']
                fsu3 = material.fsu3['fsu3']
                esu3 = material.esu3['esu3']
                esv3 = material.esv3['esv3']
                Es3 = material.Es3['Es3']
                ka3 = material.ka3['ka3']
                kb3 = material.kb3['kb3']
                kc3 = material.kc3['kc3']
                as3 = material.as3['as3']
                dm3 = material.dm3['dm3']
                psi3 = material.psi3['psi3']

                stahl4 = material.stahl4['stahl4']
                zm4 = material.zm4['zm4']
                fsy4 = material.fsy4['fsy4']
                fsu4 = material.fsu4['fsu4']
                esu4 = material.esu4['esu4']
                esv4 = material.esv4['esv4']
                Es4 = material.Es4['Es4']
                ka4 = material.ka4['ka4']
                kb4 = material.kb4['kb4']
                kc4 = material.kc4['kc4']
                as4 = material.as4['as4']
                dm4 = material.dm4['dm4']
                psi4 = material.psi4['psi4']

                Dimens=2  # Shell181
                Modell=1 #Cracked Membrane model!                                                


                # Berechnete werte
                Ec = k_E*fcc**(1/3)
                E =Ec              
                v = vc                
                
                t=t_set[count_mat]   
                nn=nn_layer_set[count_mat]                              
                h=t
                G=E/(2*(1+v))
                E111=5/6*G*t    
                E221=5/6*G*t                     
                E121=0
                
                delta_h=t/nn

                self.write_subsection(key)


            # ------------------------------------------------------------------------------------------------------
            # material properties schreiben
            # ------------------------------------------------------------------------------------------------------

            if mtype == 'ElasticIsotropic':                
                self.blank_line()
                self.write_line('mp,ex,{0},{1}'.format(count_mat+1, E))
                self.write_line('mp,ey,{0},{1}'.format(count_mat+1, E))
                self.write_line('mp,prxy,{0},{1}'.format(count_mat+1, v))
                self.write_line('mp,dens,{0},{1}'.format(count_mat+1, p))
                self.blank_line()

                self.write_line('allsel')  
                self.write_line('cmsel,s, {0}, elem'.format(ele_set[count_mat]))  
                self.write_line('mpchg,{0},all'.format(count_mat+1))

            elif mtype == 'MPCStiff':
                self.blank_line()
                self.write_line('!No Element properties for MPCs are needed') 
                self.blank_line()    

            elif mtype == 'CMMUsermat':

                self.blank_line()
                self.write_line('tb,user,{0},1,76'.format(count_mat+1))
                self.write_line('tbtemp,0')
                self.write_line('tbdata,1,{0},{1}'.format(Dimens, Modell))
                self.write_line('tbdata,3,{0},{1},{2},{3}'.format(nn, h, oo, uu))
                self.write_line('tbdata,7,{0},{1},{2},{3},{4}'.format(stahl4, zm4, as4, dm4, psi4))
                self.write_line('tbdata,12,{0},{1},{2},{3},{4}'.format(fsy4, fsu4, esu4, esv4, Es4))
                self.write_line('tbdata,17,{0},{1},{2}'.format(ka4, kb4, kc4))
                self.write_line('tbdata,20,{0},{1},{2},{3},{4}'.format(stahl3, zm3, as3, dm3, psi3))
                self.write_line('tbdata,25,{0},{1},{2},{3},{4}'.format(fsy3, fsu3, esu3, esv3, Es3))
                self.write_line('tbdata,30,{0},{1},{2}'.format(ka3, kb3, kc3))
                self.write_line('tbdata,33,{0},{1},{2},{3},{4}'.format(stahl2, zm2, as2, dm2, psi2))
                self.write_line('tbdata,38,{0},{1},{2},{3},{4}'.format(fsy2, fsu2, esu2, esv2, Es2))
                self.write_line('tbdata,43,{0},{1},{2}'.format(ka2, kb2, kc2))
                self.write_line('tbdata,46,{0},{1},{2},{3},{4}'.format(stahl1, zm1, as1, dm1, psi1))
                self.write_line('tbdata,51,{0},{1},{2},{3},{4}'.format(fsy1, fsu1, esu1, esv1, Es1))
                self.write_line('tbdata,56,{0},{1},{2}'.format(ka1, kb1, kc1))                            
                self.write_line('tbdata,59,{0},{1},{2},{3}'.format(beton, fcc, vc, ecu))            
                self.write_line('tbdata,63,{0},{1},{2},{3}'.format(k_E, theta_b0, theta_b1, k_riss))  
                self.write_line('tbdata,67,{0},{1},{2},{3},{4}'.format(lambdaTS, srmx, srmy, Begrenzung, Entfestigung))  
                self.write_line('tbdata,72,{0},{1},{2},{3}'.format(winkelD, KritQ, k_vr, fswy)) 
                self.write_line('tbdata,76,{0}'.format(R_Rohr)) 
                self.write_line('tb,state,{0},,72'.format(count_mat+1))
                self.write_line('mp,dens,{0},{1}'.format(count_mat+1, rho))
                self.blank_line()

                self.write_line('allsel')  
                self.write_line('cmsel,s, {0}, elem'.format(ele_set[count_mat]))  
                self.write_line('mpchg,{0},all'.format(count_mat+1))                   

            # ------------------------------------------------------------------------------------------------------
            # element properties schreiben
            # ------------------------------------------------------------------------------------------------------

            self.blank_line()

            # Shell
            # -------

            if mtype in ['ElasticIsotropic', 'ElasticPlastic', 'Steel', 'Concrete', 'Stiff',
                         'ConcreteSmearedCrack', 'ConcreteDamagedPlasticity', 'CMMUsermat']:

                self.write_line('sectype,{0} , shell'.format(count_mat+1))
                self.write_line('seccontrols,{0},{1},{2},,1,1,1'.format(E111, E221, E121))
                self.blank_line()
                self.write_line('*do,j,1,{0}'.format(nn))
                self.write_line('secdata,{0},{1},0,,,j'.format(delta_h,count_mat+1))
                self.write_line('secoffset,mid')
                self.write_line('*enddo')       

            # MPC
            # -------

            if mtype in ['MPCStiff']:

                self.write_line('! No Element properties for MPCs are needed') 

            # Compression
            # -----------

            if mtype in ['ElasticPlastic', 'Steel', 'Concrete', 'ConcreteSmearedCrack',
                         'ConcreteDamagedPlasticity']:

                raise NotImplementedError

            # Tension
            # -------

            if mtype in ['Concrete', 'ConcreteSmearedCrack']:

                raise NotImplementedError

            elif mtype == 'ConcreteDamagedPlasticity':

                raise NotImplementedError
            

            self.blank_line()

        self.blank_line()
        self.blank_line()











