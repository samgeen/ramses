	  @  �   k820309    �
          11.1        ��S                                                                                                           
       ../amr/bisection.f90 BISECTION                                                    
                                                          
                         @                                '�                   #NGRID    #NPART    #IGRID    #F    #U    #FP 	   #UP 
                �                                                               �                                                             �                                                                       &                                                       �                                          P                             &                   &                                                       �                                         �                 
            &                   &                                                       �                              	                                        &                   &                                                       �                             
            p                
            &                   &                                                                                                                                                            @                                                      @                                                                    &                   &                                                    @                                                   
                &                                                                                                                                              3         @                                                                    &                                                                                                            @@                                                                                          
                                                                                                                                                                                                                                              @                                                                    &                   &                                                                                                                                                                                                                                                                     
                                                      
                @                                                   
                &                   &                                                    @                                                   
                &                   &                                                    @                                                    
                &                   &                                                    @                                !                   
                &                   &                                                    @ @                               "                                   &                                                                                      #                                                       $                     @                                 %                                   &                                                                                      &                                                       '                     @                                 (                                   &                                                    @                                 )                                   &                                                    @                                 *                                   &                                                                                        +                                                       ,                                                      -                   
                &                   &                                                                                        .                                       @               64                                             /                                                                                                    0                                                       1                                   &                                                                                      2                                   &                                                    @                                 3                                   &                                                                                      4                                                       5            �                       &                                           #COMMUNICATOR                                               6            �                       &                   &                                           #COMMUNICATOR    #         @                                  7                   #X 8   #C 9   #NN :             
                                 8                   
              &                   &                                                     D                                 9                                  &                                                     
                                  :           #         @                                  ;                 #BISECTION!BUILD_BISECTION!MPI_SGI_PRIVATE_STATUS <   #BISECTION!BUILD_BISECTION!MPI_SGI_PRIVATE ?   #BISECTION!BUILD_BISECTION!MPI_SGI_PRIVATE_INPLACE B   #BISECTION!BUILD_BISECTION!MPI_SGI_PRIVATE_CHAR D   #BUILD_BISECTION%SQRT G   #BUILD_BISECTION%SUM H   #BUILD_BISECTION%MIN I   #BUILD_BISECTION%MAX J   #BUILD_BISECTION%ABS K   #BUILD_BISECTION%HUGE L   #BUILD_BISECTION%FLOOR M   #BUILD_BISECTION%DBLE N   #UPDATE O                                            <                          #BUILD_BISECTION%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE =   #BUILD_BISECTION%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE >            0�            �                  =                                 p          p            p                                           0                               >                                 p          p            p                                                                           ?                          #BUILD_BISECTION%MPI_SGI_PRIVATE%MPI_BOTTOM @   #BUILD_BISECTION%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE A            �            �                  @                             0                               A                                 p          p            p                                                                            B                          #BUILD_BISECTION%MPI_SGI_PRIVATE_INPLACE%MPI_IN_PLACE C             �            �                  C                                 p          p            p                                                                           D                          #BUILD_BISECTION%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL E   #BUILD_BISECTION%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL F   -         0�            �                  E                                 p          p            p                                  -         0                               F                                 p          p            p                                                                             G     SQRT                                            H     SUM                                            I     MIN                                            J     MAX                                            K     ABS                                            L     HUGE                                            M     FLOOR                                            N     DBLE           
                                  O           #         @                                 P                  #BISECTION!INIT_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_STATUS Q   #BISECTION!INIT_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE T   #BISECTION!INIT_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_INPLACE W   #BISECTION!INIT_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_CHAR Y   #INIT_BISECTION_HISTOGRAM%MIN \   #INIT_BISECTION_HISTOGRAM%DBLE ]                                            Q                          #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE R   #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE S            0�            �                  R                                 p          p            p                                           0                               S                                 p          p            p                                                                           T                          #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_BOTTOM U   #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE V            �            �                  U                             0                               V                                 p          p            p                                                                            W                          #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_INPLACE%MPI_IN_PLACE X             �            �                  X                                 p          p            p                                                                           Y                          #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL Z   #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL [   -         0�            �                  Z                                 p          p            p                                  -         0                               [                                 p          p            p                                                                             \     MIN                                            ]     DBLE #         @                                 ^                 #BISECTION!BUILD_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_STATUS _   #BISECTION!BUILD_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE b   #BISECTION!BUILD_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_INPLACE e   #BISECTION!BUILD_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_CHAR g   #BUILD_BISECTION_HISTOGRAM%FLOOR j   #BUILD_BISECTION_HISTOGRAM%DBLE k   #LEV l   #DIR m                                            _                          #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE `   #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE a            0�            �                  `                                 p          p            p                                           0                               a                                 p          p            p                                                                           b                          #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_BOTTOM c   #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE d            �            �                  c                             0                               d                                 p          p            p                                                                            e                          #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_INPLACE%MPI_IN_PLACE f             �            �                  f                                 p          p            p                                                                           g                          #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL h   #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL i   -         0�            �                  h                                 p          p            p                                  -         0                               i                                 p          p            p                                                                             j     FLOOR                                            k     DBLE           
                                  l                     
                                  m           %         @                                n                   
       #ROUND_TO_BISEC_RES%NINT o   #ROUND_TO_BISEC_RES%DBLE p   #X q                                              o     NINT                                            p     DBLE           
                                 q     
      #         @                                 r                  #SPLITSORT_BISECTION_HISTOGRAM%DBLE s   #LEV t   #DIR u   #WALLS v                                              s     DBLE           
                                  t                     
                                  u                     
                                 v                   
              &                                                         #BUILD_BISECTION%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE >                 #BUILD_BISECTION%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE =                 #BUILD_BISECTION%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE A                 #BUILD_BISECTION%MPI_SGI_PRIVATE%MPI_BOTTOM @                 #BUILD_BISECTION%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL F                 #BUILD_BISECTION%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL E                 #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE S                 #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE R                 #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE V                 #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_BOTTOM U                 #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL [                 #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL Z                 #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE a                 #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE `                 #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE d                 #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_BOTTOM c                 #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL i                 #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL h      �   '      fn#fn    �   @   j   AMR_PARAMETERS      @   J   AMR_COMMONS )   G  �       COMMUNICATOR+AMR_COMMONS /   �  H   a   COMMUNICATOR%NGRID+AMR_COMMONS /     H   a   COMMUNICATOR%NPART+AMR_COMMONS /   f  �   a   COMMUNICATOR%IGRID+AMR_COMMONS +   �  �   a   COMMUNICATOR%F+AMR_COMMONS +   �  �   a   COMMUNICATOR%U+AMR_COMMONS ,   R  �   a   COMMUNICATOR%FP+AMR_COMMONS ,   �  �   a   COMMUNICATOR%UP+AMR_COMMONS "   �  p       DP+AMR_PARAMETERS '     @       BISEC_ROOT+AMR_COMMONS '   Z  �       BISEC_NEXT+AMR_COMMONS '   �  �       BISEC_WALL+AMR_COMMONS $   �  q       NDIM+AMR_PARAMETERS '   �  �       BISEC_INDX+AMR_COMMONS (   �  @       NBINODES+AMR_PARAMETERS ,   �  @       NBILEAFNODES+AMR_PARAMETERS &   	  @       BOXLEN+AMR_PARAMETERS +   G	  @       ICOARSE_MAX+AMR_PARAMETERS +   �	  @       ICOARSE_MIN+AMR_PARAMETERS '   �	  @       VERBOSE+AMR_PARAMETERS !   
  @       NCPU+AMR_COMMONS '   G
  �       BISEC_HIST+AMR_COMMONS '   �
  @       BISEC_NRES+AMR_COMMONS +   +  @       NBILEVELMAX+AMR_PARAMETERS !   k  @       MYID+AMR_COMMONS &   �  @       BISEC_RES+AMR_COMMONS )   �  @       BISEC_TOL+AMR_PARAMETERS .   +  �       BISEC_CPUBOX_MIN2+AMR_COMMONS .   �  �       BISEC_CPUBOX_MAX2+AMR_COMMONS -   s  �       BISEC_CPUBOX_MIN+AMR_COMMONS -     �       BISEC_CPUBOX_MAX+AMR_COMMONS +   �  �       BISEC_CPU_LOAD+AMR_COMMONS +   G  @       JCOARSE_MIN+AMR_PARAMETERS +   �  @       KCOARSE_MIN+AMR_PARAMETERS ,   �  �       NEW_HIST_BOUNDS+AMR_COMMONS "   S  @       NX+AMR_PARAMETERS "   �  @       NY+AMR_PARAMETERS .   �  �       BISEC_HIST_BOUNDS+AMR_COMMONS +   _  �       BISEC_IND_CELL+AMR_COMMONS '   �  �       CELL_LEVEL+AMR_COMMONS $   w  @       NCOARSE+AMR_COMMONS (   �  @       NGRIDMAX+AMR_PARAMETERS    �  �       XG+AMR_COMMONS '   �  r       NVECTOR+AMR_PARAMETERS )     p       TWOTONDIM+AMR_PARAMETERS "   }  @       NZ+AMR_PARAMETERS $   �  �       CPU_MAP+AMR_COMMONS     I  �       SON+AMR_COMMONS "   �  �       FLAG1+AMR_COMMONS )   a  @       NLEVELMAX+AMR_PARAMETERS #   �  �       ACTIVE+AMR_COMMONS &   ?  �       RECEPTION+AMR_COMMONS %   �  ^       CMP_BISECTION_CPUMAP '   S  �   a   CMP_BISECTION_CPUMAP%X '   �  �   a   CMP_BISECTION_CPUMAP%C (   �  @   a   CMP_BISECTION_CPUMAP%NN     �  �      BUILD_BISECTION A   �  �   �   BISECTION!BUILD_BISECTION!MPI_SGI_PRIVATE_STATUS I   �  �      BUILD_BISECTION%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE K   &  �      BUILD_BISECTION%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE :   �  �   �   BISECTION!BUILD_BISECTION!MPI_SGI_PRIVATE ;   �  H      BUILD_BISECTION%MPI_SGI_PRIVATE%MPI_BOTTOM D   �  �      BUILD_BISECTION%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE B   o  �   �   BISECTION!BUILD_BISECTION!MPI_SGI_PRIVATE_INPLACE E   �  �      BUILD_BISECTION%MPI_SGI_PRIVATE_INPLACE%MPI_IN_PLACE ?   �   �   �   BISECTION!BUILD_BISECTION!MPI_SGI_PRIVATE_CHAR C   ^!  �      BUILD_BISECTION%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL D   "  �      BUILD_BISECTION%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL %   �"  =      BUILD_BISECTION%SQRT $   �"  <      BUILD_BISECTION%SUM $   #  <      BUILD_BISECTION%MIN $   [#  <      BUILD_BISECTION%MAX $   �#  <      BUILD_BISECTION%ABS %   �#  =      BUILD_BISECTION%HUGE &   $  >      BUILD_BISECTION%FLOOR %   N$  =      BUILD_BISECTION%DBLE '   �$  @   a   BUILD_BISECTION%UPDATE )   �$  �      INIT_BISECTION_HISTOGRAM J   L&  �   �   BISECTION!INIT_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_STATUS R   ,'  �      INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE T   �'  �      INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE C   t(  �   �   BISECTION!INIT_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE D   ?)  H      INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_BOTTOM M   �)  �      INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE K   +*  �   �   BISECTION!INIT_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_INPLACE N   �*  �      INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_INPLACE%MPI_IN_PLACE H   b+  �   �   BISECTION!INIT_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_CHAR L   5,  �      INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL M   �,  �      INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL -   }-  <      INIT_BISECTION_HISTOGRAM%MIN .   �-  =      INIT_BISECTION_HISTOGRAM%DBLE *   �-  �      BUILD_BISECTION_HISTOGRAM K   �/  �   �   BISECTION!BUILD_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_STATUS S   s0  �      BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE U   1  �      BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE D   �1  �   �   BISECTION!BUILD_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE E   �2  H      BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_BOTTOM N   �2  �      BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE L   t3  �   �   BISECTION!BUILD_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_INPLACE O   4  �      BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_INPLACE%MPI_IN_PLACE I   �4  �   �   BISECTION!BUILD_BISECTION_HISTOGRAM!MPI_SGI_PRIVATE_CHAR M   �5  �      BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL N   %6  �      BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL 0   �6  >      BUILD_BISECTION_HISTOGRAM%FLOOR /   7  =      BUILD_BISECTION_HISTOGRAM%DBLE .   D7  @   a   BUILD_BISECTION_HISTOGRAM%LEV .   �7  @   a   BUILD_BISECTION_HISTOGRAM%DIR #   �7  �       ROUND_TO_BISEC_RES (   U8  =      ROUND_TO_BISEC_RES%NINT (   �8  =      ROUND_TO_BISEC_RES%DBLE %   �8  @   a   ROUND_TO_BISEC_RES%X .   9  �       SPLITSORT_BISECTION_HISTOGRAM 3   �9  =      SPLITSORT_BISECTION_HISTOGRAM%DBLE 2   �9  @   a   SPLITSORT_BISECTION_HISTOGRAM%LEV 2   :  @   a   SPLITSORT_BISECTION_HISTOGRAM%DIR 4   Y:  �   a   SPLITSORT_BISECTION_HISTOGRAM%WALLS U   �:  P       #BUILD_BISECTION%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE-ovl#61#119 W   5;  N       #BUILD_BISECTION%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE-ovl#62#120 G   �;  I       #BUILD_BISECTION%MPI_SGI_PRIVATE%MPI_BOTTOM-ovl#64#121 P   �;  @       #BUILD_BISECTION%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE-ovl#65#122 O   <  I       #BUILD_BISECTION%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL-ovl#69#123 P   U<  H       #BUILD_BISECTION%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL-ovl#70#124 ^   �<  Y       #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE-ovl#82#125 `   �<  W       #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE-ovl#83#126 P   M=  R       #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_BOTTOM-ovl#85#127 Y   �=  I       #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE-ovl#86#128 X   �=  R       #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL-ovl#90#129 Y   :>  Q       #INIT_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL-ovl#91#130 _   �>  Z       #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUS_IGNORE-ovl#96#131 a   �>  X       #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_STATUS%MPI_STATUSES_IGNORE-ovl#97#132 Q   =?  S       #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_BOTTOM-ovl#99#133 [   �?  J       #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE%MPI_ERRCODES_IGNORE-ovl#100#134 Z   �?  S       #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGV_NULL-ovl#104#135 [   -@  R       #BUILD_BISECTION_HISTOGRAM%MPI_SGI_PRIVATE_CHAR%MPI_ARGVS_NULL-ovl#105#136 