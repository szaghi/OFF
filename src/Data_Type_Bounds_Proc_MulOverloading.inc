!> @brief Operator (*) overloading by means of type-bound procedures.
    generic:: operator(*) => self_mul_self,                                                                                        &
#ifdef r16p
                             ScalR16P_mul_self,self_mul_ScalR16P,                                                                  &
#endif
                             ScalR8P_mul_self,self_mul_ScalR8P,ScalR4P_mul_self,self_mul_ScalR4P,ScalI8P_mul_self,self_mul_ScalI8P,&
                             ScalI4P_mul_self,self_mul_ScalI4P,ScalI2P_mul_self,self_mul_ScalI2P,ScalI1P_mul_self,self_mul_ScalI1P
    ! private procedures
#ifdef r16p
    procedure, pass(self ), private:: ScalR16P_mul_self
    procedure, pass(self ), private:: self_mul_ScalR16P
#endif
    procedure, pass(self1), private:: self_mul_self
    procedure, pass(self ), private:: ScalR8P_mul_self
    procedure, pass(self ), private:: ScalR4P_mul_self
    procedure, pass(self ), private:: ScalI8P_mul_self
    procedure, pass(self ), private:: ScalI4P_mul_self
    procedure, pass(self ), private:: ScalI2P_mul_self
    procedure, pass(self ), private:: ScalI1P_mul_self
    procedure, pass(self ), private:: self_mul_ScalR8P
    procedure, pass(self ), private:: self_mul_ScalR4P
    procedure, pass(self ), private:: self_mul_ScalI8P
    procedure, pass(self ), private:: self_mul_ScalI4P
    procedure, pass(self ), private:: self_mul_ScalI2P
    procedure, pass(self ), private:: self_mul_ScalI1P
