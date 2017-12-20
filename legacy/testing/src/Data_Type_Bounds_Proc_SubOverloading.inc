!> @brief Operator (-) overloading by means of type-bound procedures.
    generic:: operator(-) => negative_self,self_sub_self,                                                                          &
#ifdef r16p
                             ScalR16P_sub_self,self_sub_ScalR16P,                                                                  &
#endif
                             ScalR8P_sub_self,self_sub_ScalR8P,ScalR4P_sub_self,self_sub_ScalR4P,ScalI8P_sub_self,self_sub_ScalI8P,&
                             ScalI4P_sub_self,self_sub_ScalI4P,ScalI2P_sub_self,self_sub_ScalI2P,ScalI1P_sub_self,self_sub_ScalI1P
    ! private procedures
#ifdef r16p
    procedure, pass(self ), private:: ScalR16P_sub_self
    procedure, pass(self ), private:: self_sub_ScalR16P
#endif
    procedure, pass(self ), private:: negative_self
    procedure, pass(self1), private:: self_sub_self
    procedure, pass(self ), private:: ScalR8P_sub_self
    procedure, pass(self ), private:: ScalR4P_sub_self
    procedure, pass(self ), private:: ScalI8P_sub_self
    procedure, pass(self ), private:: ScalI4P_sub_self
    procedure, pass(self ), private:: ScalI2P_sub_self
    procedure, pass(self ), private:: ScalI1P_sub_self
    procedure, pass(self ), private:: self_sub_ScalR8P
    procedure, pass(self ), private:: self_sub_ScalR4P
    procedure, pass(self ), private:: self_sub_ScalI8P
    procedure, pass(self ), private:: self_sub_ScalI4P
    procedure, pass(self ), private:: self_sub_ScalI2P
    procedure, pass(self ), private:: self_sub_ScalI1P
