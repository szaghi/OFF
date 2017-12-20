!> @brief Assignment overloading by means of type-bound procedures.
    generic:: assignment(=) => assign_self,    &
#ifdef r16p
                               assign_ScalR16P,&
#endif
                               assign_ScalR8P,assign_ScalR4P,assign_ScalI8P,assign_ScalI4P,assign_ScalI2P,assign_ScalI1P
    ! private procedures
#ifdef r16p
    procedure, pass(self ), private:: assign_ScalR16P
#endif
    procedure, pass(self1), private:: assign_self
    procedure, pass(self ), private:: assign_ScalR8P
    procedure, pass(self ), private:: assign_ScalR4P
    procedure, pass(self ), private:: assign_ScalI8P
    procedure, pass(self ), private:: assign_ScalI4P
    procedure, pass(self ), private:: assign_ScalI2P
    procedure, pass(self ), private:: assign_ScalI1P
