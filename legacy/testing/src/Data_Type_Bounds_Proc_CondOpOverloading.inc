!> @brief Conditional operators (=,*,/,+,-) overloading by means of type-bound procedures.
    generic:: operator(/=) => self_not_eq_self,                                                                               &
#ifdef r16p
                              R16P_not_eq_self,self_not_eq_R16P,                                                              &
#endif
                              R8P_not_eq_self,self_not_eq_R8P,R4P_not_eq_self,self_not_eq_R4P,I8P_not_eq_self,self_not_eq_I8P,&
                              I4P_not_eq_self,self_not_eq_I4P,I2P_not_eq_self,self_not_eq_I2P,I1P_not_eq_self,self_not_eq_I1P
    generic:: operator(<) => self_low_self,                                                                &
#ifdef r16p
                             R16P_low_self,self_low_R16P,                                                  &
#endif
                             R8P_low_self,self_low_R8P,R4P_low_self,self_low_R4P,I8P_low_self,self_low_I8P,&
                             I4P_low_self,self_low_I4P,I2P_low_self,self_low_I2P,I1P_low_self,self_low_I1P

    generic:: operator(<=) => self_low_eq_self,                                                                               &
#ifdef r16p
                              R16P_low_eq_self,self_low_eq_R16P,                                                              &
#endif
                              R8P_low_eq_self,self_low_eq_R8P,R4P_low_eq_self,self_low_eq_R4P,I8P_low_eq_self,self_low_eq_I8P,&
                              I4P_low_eq_self,self_low_eq_I4P,I2P_low_eq_self,self_low_eq_I2P,I1P_low_eq_self,self_low_eq_I1P
    generic:: operator(==) => self_eq_self,                                                           &
#ifdef r16p
                              R16P_eq_self,self_eq_R16P,                                              &
#endif
                              R8P_eq_self,self_eq_R8P,R4P_eq_self,self_eq_R4P,I8P_eq_self,self_eq_I8P,&
                              I4P_eq_self,self_eq_I4P,I2P_eq_self,self_eq_I2P,I1P_eq_self,self_eq_I1P
    generic:: operator(>=) => self_great_eq_self,                                                                       &
#ifdef r16p
                              R16P_great_eq_self,self_great_eq_R16P,                                                    &
#endif
                              R8P_great_eq_self,self_great_eq_R8P,R4P_great_eq_self,self_great_eq_R4P,I8P_great_eq_self,&
                              self_great_eq_I8P,I4P_great_eq_self,self_great_eq_I4P,I2P_great_eq_self,self_great_eq_I2P,&
                              I1P_great_eq_self,self_great_eq_I1P
    generic:: operator(>) => self_great_self,                                                                          &
#ifdef r16p
                             R16P_great_self,self_great_R16P,                                                          &
#endif
                             R8P_great_self,self_great_R8P,R4P_great_self,self_great_R4P,I8P_great_self,self_great_I8P,&
                             I4P_great_self,self_great_I4P,I2P_great_self,self_great_I2P,I1P_great_self,self_great_I1P
    ! private procedures
#ifdef r16p
    procedure, pass(self ), private:: R16P_not_eq_self
    procedure, pass(self ), private:: self_not_eq_R16P
    procedure, pass(self ), private:: R16P_low_self
    procedure, pass(self ), private:: self_low_R16P
    procedure, pass(self ), private:: R16P_low_eq_self
    procedure, pass(self ), private:: self_low_eq_R16P
    procedure, pass(self ), private:: R16P_eq_self
    procedure, pass(self ), private:: self_eq_R16P
    procedure, pass(self ), private:: R16P_great_eq_self
    procedure, pass(self ), private:: self_great_eq_R16P
    procedure, pass(self ), private:: R16P_great_self
    procedure, pass(self ), private:: self_great_R16P
#endif
    procedure, pass(self1), private:: self_not_eq_self
    procedure, pass(self ), private:: R8P_not_eq_self
    procedure, pass(self ), private:: R4P_not_eq_self
    procedure, pass(self ), private:: I8P_not_eq_self
    procedure, pass(self ), private:: I4P_not_eq_self
    procedure, pass(self ), private:: I2P_not_eq_self
    procedure, pass(self ), private:: I1P_not_eq_self
    procedure, pass(self ), private:: self_not_eq_R8P
    procedure, pass(self ), private:: self_not_eq_R4P
    procedure, pass(self ), private:: self_not_eq_I8P
    procedure, pass(self ), private:: self_not_eq_I4P
    procedure, pass(self ), private:: self_not_eq_I2P
    procedure, pass(self ), private:: self_not_eq_I1P
    procedure, pass(self1), private:: self_low_self
    procedure, pass(self ), private:: R8P_low_self
    procedure, pass(self ), private:: R4P_low_self
    procedure, pass(self ), private:: I8P_low_self
    procedure, pass(self ), private:: I4P_low_self
    procedure, pass(self ), private:: I2P_low_self
    procedure, pass(self ), private:: I1P_low_self
    procedure, pass(self ), private:: self_low_R8P
    procedure, pass(self ), private:: self_low_R4P
    procedure, pass(self ), private:: self_low_I8P
    procedure, pass(self ), private:: self_low_I4P
    procedure, pass(self ), private:: self_low_I2P
    procedure, pass(self ), private:: self_low_I1P
    procedure, pass(self1), private:: self_low_eq_self
    procedure, pass(self ), private:: R8P_low_eq_self
    procedure, pass(self ), private:: R4P_low_eq_self
    procedure, pass(self ), private:: I8P_low_eq_self
    procedure, pass(self ), private:: I4P_low_eq_self
    procedure, pass(self ), private:: I2P_low_eq_self
    procedure, pass(self ), private:: I1P_low_eq_self
    procedure, pass(self ), private:: self_low_eq_R8P
    procedure, pass(self ), private:: self_low_eq_R4P
    procedure, pass(self ), private:: self_low_eq_I8P
    procedure, pass(self ), private:: self_low_eq_I4P
    procedure, pass(self ), private:: self_low_eq_I2P
    procedure, pass(self ), private:: self_low_eq_I1P
    procedure, pass(self1), private:: self_eq_self
    procedure, pass(self ), private:: R8P_eq_self
    procedure, pass(self ), private:: R4P_eq_self
    procedure, pass(self ), private:: I8P_eq_self
    procedure, pass(self ), private:: I4P_eq_self
    procedure, pass(self ), private:: I2P_eq_self
    procedure, pass(self ), private:: I1P_eq_self
    procedure, pass(self ), private:: self_eq_R8P
    procedure, pass(self ), private:: self_eq_R4P
    procedure, pass(self ), private:: self_eq_I8P
    procedure, pass(self ), private:: self_eq_I4P
    procedure, pass(self ), private:: self_eq_I2P
    procedure, pass(self ), private:: self_eq_I1P
    procedure, pass(self1), private:: self_great_eq_self
    procedure, pass(self ), private:: R8P_great_eq_self
    procedure, pass(self ), private:: R4P_great_eq_self
    procedure, pass(self ), private:: I8P_great_eq_self
    procedure, pass(self ), private:: I4P_great_eq_self
    procedure, pass(self ), private:: I2P_great_eq_self
    procedure, pass(self ), private:: I1P_great_eq_self
    procedure, pass(self ), private:: self_great_eq_R8P
    procedure, pass(self ), private:: self_great_eq_R4P
    procedure, pass(self ), private:: self_great_eq_I8P
    procedure, pass(self ), private:: self_great_eq_I4P
    procedure, pass(self ), private:: self_great_eq_I2P
    procedure, pass(self ), private:: self_great_eq_I1P
    procedure, pass(self1), private:: self_great_self
    procedure, pass(self ), private:: R8P_great_self
    procedure, pass(self ), private:: R4P_great_self
    procedure, pass(self ), private:: I8P_great_self
    procedure, pass(self ), private:: I4P_great_self
    procedure, pass(self ), private:: I2P_great_self
    procedure, pass(self ), private:: I1P_great_self
    procedure, pass(self ), private:: self_great_R8P
    procedure, pass(self ), private:: self_great_R4P
    procedure, pass(self ), private:: self_great_I8P
    procedure, pass(self ), private:: self_great_I4P
    procedure, pass(self ), private:: self_great_I2P
    procedure, pass(self ), private:: self_great_I1P
