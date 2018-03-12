!> @brief Operator (+) overloading by means of type-bound procedures.
    generic:: operator(+) => positive_self,self_sum_self,                                                                          &
#ifdef r16p
                             ScalR16P_sum_self,self_sum_ScalR16P,                                                                  &
#endif
                             ScalR8P_sum_self,self_sum_ScalR8P,ScalR4P_sum_self,self_sum_ScalR4P,ScalI8P_sum_self,self_sum_ScalI8P,&
                             ScalI4P_sum_self,self_sum_ScalI4P,ScalI2P_sum_self,self_sum_ScalI2P,ScalI1P_sum_self,self_sum_ScalI1P
    ! private procedures
#ifdef r16p
    procedure, pass(self ), private:: ScalR16P_sum_self
    procedure, pass(self ), private:: self_sum_ScalR16P
#endif
    procedure, pass(self ), private:: positive_self
    procedure, pass(self1), private:: self_sum_self
    procedure, pass(self ), private:: ScalR8P_sum_self
    procedure, pass(self ), private:: ScalR4P_sum_self
    procedure, pass(self ), private:: ScalI8P_sum_self
    procedure, pass(self ), private:: ScalI4P_sum_self
    procedure, pass(self ), private:: ScalI2P_sum_self
    procedure, pass(self ), private:: ScalI1P_sum_self
    procedure, pass(self ), private:: self_sum_ScalR8P
    procedure, pass(self ), private:: self_sum_ScalR4P
    procedure, pass(self ), private:: self_sum_ScalI8P
    procedure, pass(self ), private:: self_sum_ScalI4P
    procedure, pass(self ), private:: self_sum_ScalI2P
    procedure, pass(self ), private:: self_sum_ScalI1P
