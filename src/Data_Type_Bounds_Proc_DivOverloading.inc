!> @brief Operator (/) overloading by means of type-bound procedures.
    generic:: operator(/) => self_div_self,    &
#ifdef r16p
                             self_div_ScalR16P,&
#endif
                             self_div_ScalR8P,self_div_ScalR4P,self_div_ScalI8P,self_div_ScalI4P,self_div_ScalI2P,self_div_ScalI1P
#ifdef r16p
    procedure, pass(self ), private:: self_div_ScalR16P
#endif
    procedure, pass(self1), private:: self_div_self
    procedure, pass(self ), private:: self_div_ScalR8P
    procedure, pass(self ), private:: self_div_ScalR4P
    procedure, pass(self ), private:: self_div_ScalI8P
    procedure, pass(self ), private:: self_div_ScalI4P
    procedure, pass(self ), private:: self_div_ScalI2P
    procedure, pass(self ), private:: self_div_ScalI1P
