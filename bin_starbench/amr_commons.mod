	  �3  r   k820309    �
          11.1        �'S                                                                                                           
       ../amr/amr_commons.f90 AMR_COMMONS                                                    
                                                                                                                                                                                                 d               100                                                                                                                                                        KIND                                                                                                                                                                                                                 	                                                       
                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
                                                      
                                                      
                                                      
                                                      
                                                      
                                                      
                                                      
                                                      
                                                      
                                                                                                                                                                                                                                                                                                                                             !                   
                &                                                                                     "                   
                &                                                                                     #                   
                &                                                                                     $                   
                &                                                                                        %                                                      &     
                      �                           '     d              
      p          & p        p d           p d                                          �                           (     d              
      p          & p        p d           p d                                          �                           )     d              
      p          & p        p d           p d                                          �                           *     d              
      p          & p        p d           p d                                          �                           +     d              
      p          & p        p d           p d                                          �                           ,     d              
      p          & p        p d           p d                                          �                           -     d              
      p          & p        p d           p d                                          �                            .     d                    p          & p        p d           p d                                          �                            /     d                    p          & p        p d           p d                                          �                            0     d                    p          & p        p d           p d                                          �                           1     d              
      p          & p        p d           p d                                          �                           2     d              
      p          & p        p d           p d                                          �                           3     d              
      p          & p        p d           p d                                        �                            4     d                    p          & p        p d           p d                                                                    5                                   &                   &                                                                                      6                                   &                   &                                                                                      7                                   &                   &                                                                                      8                                   &                   &                                                                                      9                                   &                   &                                                                                      :                                   &                   &                                                                                      ;                                   &                   &                                                                                        <                                                         =                                                         >                                                         ?                                                         @                                                      A                   
                &                   &                                                                                      B                                   &                   &                                                                                      C                                   &                                                                                      D                                   &                                                                                      E                                   &                                                                                      F                                   &                                                                                      G                                   &                                                                                      H                                   &                                                                                      I                                   &                                                                                      J                                   &                                                                                     K                   
                &                                                                                     L                   
                &                                                                                     M                   
                &                                                                                       N     
                                                   O     
                                                 P                   
                &                                                                                      Q                                   &                   &                                                                                        R                                                       S                                   &                                                                                     T                   
                &                   &                                                                                     U                   
                &                   &                                                                                     V                   
                &                   &                                                                                     W                   
                &                   &                                                                                      X                                   &                                                                                      Y                                   &                   &                                                                                      Z                                   &                                                                                      [                                   &                                                                                      \                                   &                                                                                      ]                                   &                                                                                       ^     
                                                    _                              @                           `     '�                   #NGRID a   #NPART b   #IGRID c   #F d   #U e   #FP f   #UP g                �                               a                                �                               b                              �                              c                                         &                                                       �                              d            P                             &                   &                                                       �                             e            �                 
            &                   &                                                       �                              f                                        &                   &                                                       �                             g            p                
            &                   &                                                                                      h            �                       &                                           #COMMUNICATOR `                                              i            �                       &                   &                                           #COMMUNICATOR `                                              j            �                       &                   &                                           #COMMUNICATOR `                                              k            �                       &                   &                                           #COMMUNICATOR `                                             l                                                      m                                                      n                                                      o     
                                                 p     
                                                 q     
          �   +      fn#fn    �   @   J   AMR_PARAMETERS "     p       DP+AMR_PARAMETERS (   {  s       MAXLEVEL+AMR_PARAMETERS #   �  p       QDP+AMR_PARAMETERS $   ^  =       KIND+AMR_PARAMETERS    �  @       OUTPUT_DONE    �  @       INIT      @       BALANCE    [  @       SHRINK    �  @       NSTEP    �  @       NSTEP_COARSE !     @       NSTEP_COARSE_OLD    [  @       NFLAG    �  @       NCREATE    �  @       NKILL      @       NCOARSE    [  @       NGRID_CURRENT    �  @       EMAG_TOT    �  @       EKIN_TOT      @       EINT_TOT    [  @       EPOT_TOT    �  @       EPOT_TOT_OLD    �  @       EPOT_TOT_INT      @       CONST    [  @       AEXP_OLD    �  @       RHO_TOT    �  @       T      @       NCPU    [  @       NDOMAIN    �  @       MYID    �  @       OVERLOAD    	  @       N_FRW    [	  �       AEXP_FRW    �	  �       HEXP_FRW    s
  �       TAU_FRW    �
  �       T_FRW    �  @       NLEVELMAX_PART    �  @       AEXP_INI      �       DFACT    �  �       ASTART    S  �       VFACT    �  �       XOFF1    �  �       XOFF2    ?  �       XOFF3    �  �       DXINI    �  �       N1    +  �       N2    �  �       N3    s  �       DTOLD      �       DTNEW    �  �       RHO_MAX    _  �       NSUBCYCLE      �       HEADL    �  �       TAILL    K  �       NUMBL    �  �       NUMBTOT    �  �       HEADB    7  �       TAILB    �  �       NUMBB      @       HEADF    �  @       TAILF    �  @       NUMBF    ?  @       USED_MEM      @       USED_MEM_TOT    �  �       XG    c  �       NBOR      �       FATHER    �  �       NEXT      �       PREV    �  �       SON    7  �       FLAG1    �  �       FLAG2    O  �       CPU_MAP    �  �       CPU_MAP2    g   �       HILBERT_KEY    �   �       BOUND_KEY    !  �       BOUND_KEY2    "  @       ORDER_ALL_MIN    K"  @       ORDER_ALL_MAX    �"  �       BISEC_WALL    #  �       BISEC_NEXT    �#  @       BISEC_ROOT    �#  �       BISEC_INDX !   �$  �       BISEC_CPUBOX_MIN !   +%  �       BISEC_CPUBOX_MAX "   �%  �       BISEC_CPUBOX_MIN2 "   s&  �       BISEC_CPUBOX_MAX2    '  �       BISEC_CPU_LOAD    �'  �       BISEC_HIST "   G(  �       BISEC_HIST_BOUNDS     �(  �       NEW_HIST_BOUNDS    _)  �       BISEC_IND_CELL    �)  �       CELL_LEVEL    w*  @       BISEC_RES    �*  @       BISEC_NRES    �*  �       COMMUNICATOR #   �+  H   a   COMMUNICATOR%NGRID #   �+  H   a   COMMUNICATOR%NPART #   ,  �   a   COMMUNICATOR%IGRID    �,  �   a   COMMUNICATOR%F    V-  �   a   COMMUNICATOR%U     .  �   a   COMMUNICATOR%FP     �.  �   a   COMMUNICATOR%UP    Z/  �       ACTIVE    �/  �       BOUNDARY    �0  �       EMISSION    d1  �       RECEPTION    2  @       TYPE_HYDRO    Z2  @       TYPE_ACCEL    �2  @       TYPE_FLAG    �2  @       UNITS_DENSITY    3  @       UNITS_TIME    Z3  @       UNITS_LENGTH 