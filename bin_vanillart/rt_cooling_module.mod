	  �  :   k820309    �
          11.1        NK(S                                                                                                           
       ../rt/rt_cooling_module.f90 RT_COOLING_MODULE              RT_SET_MODEL RT_SOLVE_COOLING UPDATE_UVRATES CMP_CHEM_EQ ISHE X Y RHOC KB MH T2_MIN_FIX TWOPI N_U INPU IFPU SIGNC SIGEC PHRATE UVRATES                                                     
       MYID                                                     
       X Y                                                    
                                                          
                      !                                
                      !                                
       #         @                                                     #NMODEL    #J0IN_IN 	   #J0MIN_IN 
   #ALPHA_IN    #NORMFACJ0_IN    #ZREIONIZ_IN    #CORRECT_COOLING    #REALISTIC_NE    #H    #OMEGAB    #OMEGA0    #OMEGAL    #ASTART_SIM    #T2_SIM                @                                                      @                             	     
                   @                             
     
                   @                                  
                   @                                  
                   @                                  
                   @                                                      @                                                    D @@                                  
                 D @@                                  
                 D @@                                  
                 D @@                                  
                   @                                  
                 D  @                                  
       #         @                                                	   #RT_SOLVE_COOLING%MIN    #RT_SOLVE_COOLING%MAX    #U    #DNPDT    #DFPDT    #NH    #C_SWITCH    #ZSOLAR    #DT    #A_EXP     #NCELL !                 @                                 MIN               @                                 MAX           D @@  �                                �             
 	    p �        & p        p �        p            p �        p                                    D @@  �                                �             
 
    p �        & p        p �        p            p �        p                                    D @@  �                                �             
     p �        & p        p �        p            p �        p                                    D @@  �                                �             
     p          & p        p �          p �                                  D @@  �                                 �                  p          & p        p �          p �                                  D @@  �                                �             
     p          & p        p �          p �                                  D @@                                  
                 D @@                                   
                   @                              !            #         @                                  "                    #         @                                  #                  #CMP_CHEM_EQ%ABS $   #CMP_CHEM_EQ%MAX %   #TK &   #NH '   #T_RAD_SPEC (   #NSPEC )   #NTOT *   #MU +                 @                            $     ABS               @                            %     MAX           
  @@                             &     
                
   @                             '     
                
   @  �                           (                   
 "   p          & p        p            p                                    D  @  �                           )                   
 #    p          & p        p            p                                    D  @                             *     
                 D  @                             +     
                 @                                ,                                                        -     
                 
                 {�!����9        1.88000D-29                                            .     
                 
                 �����<        1.38062D-16                                            /     
                 
                 �W��� ;        1.66000D-24                                            0     
                 
                 {�G�z�?        1.D-2                                            1     
                 
                 ���S�!@        6.2831853D0                                             2                                                                    @                                 3                         p          p            p                                     @                                 4                         p          p            p                                                                      5                   
      p          p          p            p          p                                                                      6                   
      p          p          p            p          p                                                                      7                   
      p          p          p            p          p                                     @ @                              8                   
      p          p          p            p          p                             �   6      fn#fn '   �   �   b   uapp(RT_COOLING_MODULE    m  E   J  AMR_COMMONS    �  D   J  COOLING_MODULE    �  @   J   RT_PARAMETERS !   6  @   J   COOLRATES_MODULE !   v  @       X+COOLING_MODULE !   �  @       Y+COOLING_MODULE    �        RT_SET_MODEL $     @   a   RT_SET_MODEL%NMODEL %   D  @   a   RT_SET_MODEL%J0IN_IN &   �  @   a   RT_SET_MODEL%J0MIN_IN &   �  @   a   RT_SET_MODEL%ALPHA_IN *     @   a   RT_SET_MODEL%NORMFACJ0_IN )   D  @   a   RT_SET_MODEL%ZREIONIZ_IN -   �  @   a   RT_SET_MODEL%CORRECT_COOLING *   �  @   a   RT_SET_MODEL%REALISTIC_NE      @   a   RT_SET_MODEL%H $   D  @   a   RT_SET_MODEL%OMEGAB $   �  @   a   RT_SET_MODEL%OMEGA0 $   �  @   a   RT_SET_MODEL%OMEGAL (     @   a   RT_SET_MODEL%ASTART_SIM $   D  @   a   RT_SET_MODEL%T2_SIM !   �  �       RT_SOLVE_COOLING %   ]  <      RT_SOLVE_COOLING%MIN %   �  <      RT_SOLVE_COOLING%MAX #   �  �   a   RT_SOLVE_COOLING%U '   �	  �   a   RT_SOLVE_COOLING%DNPDT '   ]
  �   a   RT_SOLVE_COOLING%DFPDT $   !  �   a   RT_SOLVE_COOLING%NH *   �  �   a   RT_SOLVE_COOLING%C_SWITCH (   i  �   a   RT_SOLVE_COOLING%ZSOLAR $     @   a   RT_SOLVE_COOLING%DT '   M  @   a   RT_SOLVE_COOLING%A_EXP '   �  @   a   RT_SOLVE_COOLING%NCELL    �  H       UPDATE_UVRATES      �       CMP_CHEM_EQ     �  <      CMP_CHEM_EQ%ABS        <      CMP_CHEM_EQ%MAX    <  @   a   CMP_CHEM_EQ%TK    |  @   a   CMP_CHEM_EQ%NH '   �  �   a   CMP_CHEM_EQ%T_RAD_SPEC "   `  �   a   CMP_CHEM_EQ%NSPEC !     @   a   CMP_CHEM_EQ%NTOT    D  @   a   CMP_CHEM_EQ%MU    �  @       ISHE    �  {       RHOC    ?  {       KB    �  {       MH    5  u       T2_MIN_FIX    �  {       TWOPI    %  p       N_U    �  �       INPU    )  �       IFPU    �  �       SIGNC    q  �       SIGEC    %  �       PHRATE    �  �       UVRATES 