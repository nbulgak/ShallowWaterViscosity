#define _CRT_SECURE_NO_WARNINGS
#include"ShallowWaterViscosity.h"
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ Godunov block*/
static double dfdp(double P, double pressure, double density)
{
	double Pi, SqrtOfPi, c;

	Pi = P / pressure;
	c = sqrt(gamma*(pressure / density)); // value of sonic velocity
	SqrtOfPi = sqrt(((gamma + 1.) / (2.*gamma))*Pi + ((gamma - 1.) / (2.*gamma)));

	if (P < pressure)
		return (1. / (gamma*P)) * c*pow(Pi, ((gamma - 1.) / (2.*gamma)));
	else
		return ((gamma + 1.)*Pi + (3.*gamma - 1.)) / (4.*gamma*density*c * SqrtOfPi*SqrtOfPi*SqrtOfPi);
}
static double f(double P, double pressure, double density)
{
	double Pi, c;

	Pi = P / pressure;
	c = sqrt(gamma*(pressure / density)); // value of sonic velocity

	if (P < pressure)
		return (2.*c / (gamma - 1.))*(pow(Pi, ((gamma - 1.) / (2.*gamma))) - 1.);
	else
		return (P - pressure) / (density*c * sqrt(((gamma + 1.) / (2.*gamma))*Pi + ((gamma - 1.) / (2.*gamma))));
}
static double F1(double P, double pressureI, double densityI, double pressureII, double densityII)
{
	return f(P, pressureI, densityI) + f(P, pressureII, densityII);
}
static double RootSeeker(double P0, double eps, struct Values* ValI, struct Values* ValII, double a, double b)
{
	double Pnext = -1., P, k, Pprev;
	double r, c, d, fc, fd;
	r = (double)(3. - sqrt(5.)) / 2.;
	P = P0;

	if (P>0)
		Pnext = P - (F1(P, ValI->p, ValI->rho, ValII->p, ValII->rho) - (ValI->u - ValII->u)) / (dfdp(P, ValI->p, ValI->rho) + dfdp(P, ValII->p, ValII->rho));

	while ((fabs(Pnext - P) > eps) && (!(Pnext < 0)))
	{
		Pprev = P;
		P = Pnext;
		Pnext = P - (F1(P, ValI->p, ValI->rho, ValII->p, ValII->rho) - (ValI->u - ValII->u)) / (dfdp(P, ValI->p, ValI->rho) + dfdp(P, ValII->p, ValII->rho));
		// printf("Pnext = %.15lf\n",Pnext);
		if (((P<Pprev) && (Pprev<Pnext)) || ((P>Pprev) && (Pprev>Pnext)))
			Pnext = (Pnext + P) / 2;
		if (!(fabs(Pnext - Pprev) > eps))
			Pnext = (Pnext + P) / 2;
	}

	if ((Pnext < 0))
	{
		c = a + r*(b - a);
		d = b - r*(b - a);

		fc = F1(c, ValI->p, ValI->rho, ValII->p, ValII->rho) - (ValI->u - ValII->u);
		fd = F1(d, ValI->p, ValI->rho, ValII->p, ValII->rho) - (ValI->u - ValII->u);

		while (b - a > eps*0.5)
		{
			k = b - a;
			if (fd > 0)
			{
				b = d;
				d = c;
				c = a + r*(b - a);
				fd = fc;
				fc = F1(c, ValI->p, ValI->rho, ValII->p, ValII->rho) - (ValI->u - ValII->u);
				//              printf("Goldy: fc\n");
			}
			if (fc < 0)
			{
				a = c;
				c = d;
				d = b - r*(b - a);
				fc = fd;
				fd = F1(d, ValI->p, ValI->rho, ValII->p, ValII->rho) - (ValI->u - ValII->u);
				//            printf("Goldy: fd\n");
			}
			//   printf ("k=%.15lf\n",k);
		}

		P = (b + a) / 2.;

	}

	return P;
}
static double QuasyMassVelocity(double P, struct Values *Val)
{
	double Pi, c;

	Pi = P / Val->p;//pressure ratio
	c = sqrt(gamma*(Val->p / Val->rho)); // value of sonic velocity

	if (P < Val->p)
		return ((gamma - 1.) / (2.*gamma)) * Val->rho * c * ((1. - Pi) / (1. - pow(Pi, ((gamma - 1.) / (2.*gamma)))));
	else
		return sqrt(Val->rho*(P*((gamma + 1.) / 2.) + Val->p*((gamma - 1.) / 2.)));
}
static void Crumbling(struct Values*v1, struct Values*v2, struct Values*bv, double VelocityX)
{
	double eps = 0.000000001;
	double Reverse, UpIIepsU, U10pII, U11pII, UpIIepsD, UpIepsU, UpIepsD, Up0epsU, Up0, UpI,
		P,
		a1, aI, a2, aII,
		cI, cII, cIStar, cIIStar,
		TettaI, TettaII,
		UOfContact, DI, DII, DIStar, DIIStar, cStar;
	int IndicOfLinear = 0, IndicOfAssignment = 0, IndicOfRevers = 0;

	struct Values *ValI, *ValII, *BigVal;
	ValI = (struct Values*) malloc(sizeof(struct Values));
	ValII = (struct Values*) malloc(sizeof(struct Values));
	BigVal = (struct Values*) malloc(sizeof(struct Values));

	ValI->u = v1->u;
	ValII->u = v2->u;

	ValI->p = v1->p;
	ValII->p = v2->p;
	ValI->rho = v1->rho;
	ValII->rho = v2->rho;

	UOfContact = ValI->u;


	//    printf("P1 = %lf,   P2 = %lf\n", ValI->p, ValII->p);
	if ((ValI->p) > (ValII->p))
	{
		//        printf("Revers\n");
		Reverse = ValI->p;
		ValI->p = ValII->p;
		ValII->p = Reverse;

		Reverse = ValI->rho;
		ValI->rho = ValII->rho;
		ValII->rho = Reverse;

		Reverse = ValI->u;
		ValI->u = -ValII->u;
		ValII->u = -Reverse;

		VelocityX = -VelocityX;

		IndicOfRevers++;
		//        printf("Indicator of revers = %d\n", IndicOfRevers);

	}

	UpIIepsU = F1(ValII->p + eps, ValI->p, ValI->rho, ValII->p, ValII->rho);//skorost' razleta + eps(?????????????????)

	if (ValI->u - ValII->u > UpIIepsU)//Ударные волны
	{
		//        printf("Shock waves\n");

		U10pII = F1(ValII->p * 10., ValI->p, ValI->rho, ValII->p, ValII->rho);//skorost' hui znaet chego(linearizacii)

		if ((ValI->u - ValII->u) > U10pII)//Считаем, что P будет много больше pI и pII
		{
			TettaI = (gamma + 1.) / 2.;
			TettaII = TettaI;

			a1 = 1 / sqrt(TettaI * ValI->rho) + 1 / sqrt(TettaII * ValII->rho);
			a2 = ValI->p / sqrt(TettaI * ValI->rho) + ValII->p / sqrt(TettaII * ValII->rho);

			P = (ValI->u - ValII->u + sqrt((ValI->u - ValII->u)*(ValI->u - ValII->u) + 4 * a1*a2)) / (2 * a1);
			P = P*P;

		}
		else
		{
			U11pII = F1(ValII->p * 1.1, ValI->p, ValI->rho, ValII->p, ValII->rho);// ONA POCHTI P2 (linearizacii)

			if (ValI->u - ValII->u > U11pII)//P сравнимо с pII
			{
				cI = sqrt(gamma*(ValI->p / ValI->rho)); // value of sonic velocity
				cII = sqrt(gamma*(ValII->p / ValII->rho)); // value of sonic velocity

				P = (ValI->p * ValII->rho*cII + ValII->p * ValI->rho*cI + (ValI->u - ValII->u)*ValI->rho*cI*ValII->rho*cII) / (ValII->rho*cII + ValI->rho*cI);

			}
			else
			{
				if (ValII->p > 10 * ValI->p)
				{
					TettaII = gamma;
					TettaI = (gamma + 1) / 2;
					IndicOfLinear = 1;

				}
				else
				{
					if (ValII->p > 1.1*ValI->p)
					{
						cI = sqrt(gamma*(ValI->p / ValI->rho)); // value of sonic velocity
						cII = sqrt(gamma*(ValII->p / ValII->rho)); // value of sonic velocity

						P = (ValI->p * ValII->rho*cII + ValII->p * ValI->rho*cI + (ValI->u - ValII->u)*ValI->rho*cI*ValII->rho*cII) / (ValII->rho*cII + ValI->rho*cI);

					}
					else
					{
						TettaI = gamma;
						TettaII = TettaI;
						IndicOfLinear = 1;
					}
				}
			}

			if (IndicOfLinear == 1)
			{
				a1 = 1 / sqrt(TettaI * ValI->rho) + 1 / sqrt(TettaII * ValII->rho);
				a2 = ValI->p / sqrt(TettaI * ValI->rho) + ValII->p / sqrt(TettaII * ValII->rho);

				P = (ValI->u - ValII->u + sqrt((ValI->u - ValII->u)*(ValI->u - ValII->u) + 4 * a1*a2)) / (2 * a1);
				P = P*P;
			}
			if (F1(P, ValI->p, ValI->rho, ValII->p, ValII->rho) > (ValI->u - ValII->u))
			{
				if (F1(P - eps, ValI->p, ValI->rho, ValII->p, ValII->rho) > (ValI->u - ValII->u))
					P = RootSeeker(P, eps, ValI, ValII, ValII->p + eps, P);
			}
			else
			{
				if (F1(P + eps, ValI->p, ValI->rho, ValII->p, ValII->rho) < (ValI->u - ValII->u))
					P = RootSeeker(P, eps, ValI, ValII, P, 10 * ValII->p);
			}
		}
		//    printf("Shock waves\n");
	}
	else
	{
		UpIIepsD = F1(ValII->p - eps, ValI->p, ValI->rho, ValII->p, ValII->rho);

		if (ValI->u - ValII->u > UpIIepsD)// Удар и однородное состояние
		{
			P = ValII->p;
		}
		else
		{
			UpIepsU = F1(ValI->p + eps, ValI->p, ValI->rho, ValII->p, ValII->rho);

			if (ValI->u - ValII->u > UpIepsU)//Shock wave and expansion fan
			{
				//                printf("Shock wave and expansion fan\n");

				cI = sqrt(gamma*(ValI->p / ValI->rho)); // value of sonic velocity
				cII = sqrt(gamma*(ValII->p / ValII->rho)); // value of sonic velocity

				P = (ValI->p * ValII->rho*cII + ValII->p * ValI->rho*cI + (ValI->u - ValII->u)*ValI->rho*cI*ValII->rho*cII) / (ValII->rho*cII + ValI->rho*cI);

				P = RootSeeker(P, eps, ValI, ValII, ValI->p + eps, ValII->p - eps);

			}
			else
			{
				UpIepsD = F1(ValI->p - eps, ValI->p, ValI->rho, ValII->p, ValII->rho);

				if (ValI->u - ValII->u >UpIepsD)// Волна разрежения и однородное состояние
				{
					P = ValI->p;
				}
				else
				{
					Up0epsU = F1(eps, ValI->p, ValI->rho, ValII->p, ValII->rho); //U of p0+eps
					//                    printf("Up0epsU = %lf\n", Up0epsU);
					if (ValI->u - ValII->u > Up0epsU)//Two expansion fans
					{
						//                        printf("Two expansion fans\n");

						Up0 = F1(0, ValI->p, ValI->rho, ValII->p, ValII->rho);
						UpI = F1(ValI->p, ValI->p, ValI->rho, ValII->p, ValII->rho);

						P = ValI->p * pow(((ValI->u - ValII->u - Up0) / (UpI - Up0)), 2.*gamma / (gamma - 1.));
					}
					else //Vacuum
					{
						P = 0.;
						//                        printf("Vacuum\n");
					}
				}
			}
		}
	}
	//    printf("P=%lf\n",P);
	/*    Reverse = Test(P,ValI,ValII);
	printf("F(P)=%lf\n",Reverse);*/

	cI = sqrt(gamma*(ValI->p / ValI->rho)); // value of sonic velocity
	cII = sqrt(gamma*(ValII->p / ValII->rho)); // value of sonic velocity
	/*
	//Time searching
	if((P > ValI->p)&&(P > ValII->p))
	{
	aI = QuasyMassVelocity( P, ValI);
	aII = QuasyMassVelocity( P, ValII);

	DII = ValII->u + aII/ValII->rho;
	DI = ValI->u - aI/ValI->rho;
	//        printf("1\n");
	}
	if((P > ValI->p)&&(P < ValII->p))
	{
	aI = QuasyMassVelocity( P, ValI);

	DII = ValII->u + cII;
	DI = ValI->u - aI/ValI->rho;
	//        printf("2\n");
	}
	if((P < ValI->p)&&(P < ValII->p))
	{
	//        printf("3\n");

	if((ValI->p > 0)&&(ValII->p > 0))
	{
	DII = ValII->u + cII;
	DI = ValI->u - cI;
	}
	else
	{
	if(ValI->p > 0)
	{
	DII = ValI->u + 2.*cI/(gamma - 1.);
	DI = ValI->u - cI;
	}
	if(ValII->p > 0)
	{
	DII = ValII->u + cII;
	DI = ValII->u - 2.*cII/(gamma - 1.);
	}
	}
	}
	//    printf("DI = %lf    DII = %lf\n", DI, DII);

	if(!( fabs(h/DI) < fabs(h/DII) ))
	Tau = fabs(h/DII);
	else
	Tau = fabs(h/DI);
	//Time searching
	*/
	while (IndicOfAssignment == 0)//Search big values
	{
		if (!(P < ValII->p))
		{
			aII = QuasyMassVelocity(P, ValII);
			DII = ValII->u + aII / ValII->rho;

			//            printf("DII=%lf\n",DII);
			aI = QuasyMassVelocity(P, ValI);
			UOfContact = (aI*ValI->u + aII*ValII->u + ValI->p - ValII->p) / (aI + aII);

			if (!(VelocityX < DII))
			{

				//                printf("Area of initial configuration\n");

				BigVal->p = ValII->p;
				BigVal->u = ValII->u;
				BigVal->rho = ValII->rho;

				//   UOfContact = (aI*ValI->u + aII*ValII->u + ValI->p - ValII->p)/(aI+aII);

				IndicOfAssignment = 1;
			}
			else
			{
				//  aI = QuasyMassVelocity( P, ValI);
				// UOfContact = (aI*ValI->u + aII*ValII->u + ValI->p - ValII->p)/(aI+aII);
				//                printf("Velocity of contact = %lf\n", UOfContact);

				if (!(VelocityX < UOfContact))
				{

					//                    printf("Area of gas diturbed by shock wave\n");

					BigVal->p = P;
					BigVal->u = UOfContact;
					BigVal->rho = ValII->rho*aII / (aII + ValII->rho*(ValII->u - UOfContact));
					IndicOfAssignment = 1;
				}
			}
		}
		else//2 exp fan
		{
			DII = ValII->u + cII;
			aI = QuasyMassVelocity(P, ValI);
			aII = QuasyMassVelocity(P, ValII);
			UOfContact = (aI*ValI->u + aII*ValII->u + ValI->p - ValII->p) / (aI + aII);

			if (!(VelocityX < DII))
			{

				//                printf("Area of initial configuration\n");

				BigVal->p = ValII->p;
				BigVal->u = ValII->u;
				BigVal->rho = ValII->rho;

				//  UOfContact = (aI*ValI->u + aII*ValII->u + ValI->p - ValII->p)/(aI+aII);

				IndicOfAssignment = 1;
			}
			else
			{
				if (P > 0)
				{
					//  aI = QuasyMassVelocity( P, ValI);
					// aII = QuasyMassVelocity( P, ValII);

					//   UOfContact = (aI*ValI->u + aII*ValII->u + ValI->p - ValII->p)/(aI+aII);

					cIIStar = cII - ((gamma - 1.) / 2.)*(ValII->u - UOfContact);
					DIIStar = UOfContact + cIIStar;

					if (!(VelocityX < DIIStar))
					{

						//                        printf("Area of expansion fan\n");

						cStar = -((gamma - 1.) / (gamma + 1.))*(ValII->u - VelocityX) + 2.*cII / (gamma + 1.);
						//                        printf("cStar = %lf\n", cStar);

						BigVal->u = VelocityX - cStar;
						BigVal->p = ValII->p * pow((cStar / cII), (2.*gamma / (gamma - 1.)));
						BigVal->rho = BigVal->p * gamma / (cStar*cStar);
						IndicOfAssignment = 1;
					}
					else
					{
						if (!(VelocityX < UOfContact))
						{
							//                            printf("Area of gas disturbed by expansion fan\n");

							BigVal->p = P;
							BigVal->u = UOfContact;
							BigVal->rho = gamma*(P / (cIIStar*cIIStar));
							IndicOfAssignment = 1;
						}
					}
				}
				else
				{

					if (!(VelocityX < ValII->u - 2.*cII / (gamma - 1.)))
					{
						//                        printf("Area of expansion fan at vacuum\n");

						cStar = -((gamma - 1.) / (gamma + 1.))*(ValII->u - VelocityX) + 2.*cII / (gamma + 1.);

						BigVal->p = ValII->p*pow(cStar / cII, (2.*gamma / (gamma - 1.)));
						BigVal->rho = gamma*BigVal->p / (cStar*cStar);
						BigVal->u = VelocityX - cStar;
						IndicOfAssignment = 1;
					}
					else
					{
						if (!(VelocityX > ValI->u + 2.*cI / (gamma - 1.)))
							//                        printf("Area of vacuum\n");

							BigVal->p = 0;
						BigVal->u = 0;
						BigVal->rho = 0;
						IndicOfAssignment = 1;
					}
				}
			}
		}

		if (IndicOfAssignment == 0)
		{

			//            printf("Assignment\n");

			Reverse = ValI->p;
			ValI->p = ValII->p;
			ValII->p = Reverse;

			Reverse = ValI->rho;
			ValI->rho = ValII->rho;
			ValII->rho = Reverse;

			Reverse = ValI->u;
			ValI->u = -ValII->u;
			ValII->u = -Reverse;

			VelocityX = -VelocityX;
			Reverse = cI;
			cI = cII;
			cII = Reverse;

			//            printf("Indicator of revers = %d\n", IndicOfRevers);

			IndicOfRevers = (IndicOfRevers + 1) % 2;
			//            printf("Indicator of revers = %d\n", IndicOfRevers);

		}
	}

	if (IndicOfRevers == 1)
	{
		/*        Reverse = ValI->p;
		ValI->p = ValII->p;
		ValII->p = Reverse;

		Reverse = ValI->rho;
		ValI->rho = ValII->rho;
		ValII->rho = Reverse;

		Reverse = ValI->u;
		ValI->u = -ValII->u;
		ValII->u = -Reverse;*/

		VelocityX = -VelocityX;

		BigVal->u = -BigVal->u;
		UOfContact = -UOfContact;
	}

	bv->u = BigVal->u;
	bv->rho = BigVal->rho;
	bv->p = BigVal->p;


	free(ValI); free(ValII); free(BigVal);
	return;
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void Create_Square(All_data *data, int a, int b)
{
	int i, j;

	for (i = (a - 1)*data->n + a; i <= (a - 1)*data->n + b - 1; ++i)
	{
		data->layer[0][i].c = 'u';
		data->layer[1][i].c = 'u';
	}

	for (i = b*data->n + a; i <= b*data->n + b - 1; ++i)
	{
		data->layer[0][i].c = 'd';
		data->layer[1][i].c = 'd';
	}
	
	for (i = a*data->n + a - 1; i < b*data->n + a - 1; i+=data->n )
	{
		data->layer[0][i].c = 'r';
		data->layer[1][i].c = 'r';
	}

	for (i = a*data->n + b; i < b*data->n + b; i += data->n)
	{
		data->layer[0][i].c = 'l';
		data->layer[1][i].c = 'l';
	}

	for (i = a; i < b; ++i)
	{
		for (j = a; j < b; ++j)
		{
			data->layer[1][j + i*data->n].square = 'y';
			data->layer[0][j + i*data->n].square = 'y';
		}
	}
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void B2F(struct Values1* val, double* F, double Fr)
{
	F[0] = val->h*val->u;
	F[1] = val->h*val->u*val->u + val->h*val->h / Fr / Fr / 2.;/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
	F[2] = val->h*val->u*val->v;
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void B2E(struct Values1* val, double* E, double Fr)
{
	E[0] = val->h*val->v;
	E[1] = val->h*val->v*val->u;
	E[2] = val->h*val->v*val->v + val->h*val->h / Fr / Fr / 2.;/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static int Point_check(Point* p)
{
	if (p->h < MACHEPS)
		return 1;
	if (p->h != p->h)
		return 2;
	if (p->u != p->u)
		return 2;
	if (p->v != p->v)
		return 2;
	return 0;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void u2U(Point* p)
{
	p->u1 = p->h;
	p->u2 = p->h*p->u;
	p->u3 = p->h*p->v;
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void U2u(Point* p)
{

	p->h = p->u1;
	if (p->h > MACHEPS)
	{
		p->u = p->u2 / p->h;
		p->v = p->u3 / p->h;
	}
	else
		printf("bad h!\n");
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
/*p2->p1*/
static void assign_points(Point* p1, Point* p2) 
{
	p1->h = p2->h; p1->u = p2->u; p1->v = p2->v;
	p1->u1 = p2->u1; p1->u2 = p2->u2; p1->u3 = p2->u3;
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void Print_grid(All_data *data)
{
	FILE* out = fopen("grid.dat", "w");
	int i, j;
	unsigned int n = data->n, m = data->m;
	fprintf(out, "FILETYPE = GRID\n");
	fprintf(out, "VARIABLES = \"x\" \"y\" \n");
	fprintf(out, "ZONE\nI=%d, J=%d\n ZONETYPE = Ordered, DATAPACKING=BLOCK\n", n+1, m+1);

	for (i = 0; i < m + 1 ; ++i)//output x
	{
		for (j = 0; j < n; ++j)
		{
			if (i == m)
			{
				fprintf(out, "%lf ", data->layer[1][j + (i - 1)*n].x);
			}
			else
			    fprintf(out, "%lf ", data->layer[1][j + i*n].x);
		}
		if (i == m)
			fprintf(out, "%lf\n", (data->layer[1][j + (i - 1)*n - 1].x + data->d_x));
		else
		    fprintf(out, "%lf\n", (data->layer[1][j + i*n - 1].x + data->d_x));
	}
	//fprintf(out, "\n");
	for (i = 0; i < m + 1; ++i)//output y
	{
		for (j = 0; j < n; ++j)
		{
			if (i == m)
				fprintf(out, "%lf ", data->layer[1][j + (i - 1)*n].y + data->d_y);
			else
			    fprintf(out, "%lf ", data->layer[1][j + i*n].y);
		}
		if (i == m)
			fprintf(out, "%lf\n", data->layer[1][j + (i - 1)*n - 1].y + data->d_y);
		else
		    fprintf(out, "%lf\n", (data->layer[1][j + i*n - 1].y));
	}

	fclose(out);
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void Print_plotina(All_data *data)
{
	char s[64];
	int i, j, k;
	double light_vel = sqrt(10. * 9.80665);
	unsigned int n = data->n, m = data->m;
	sprintf(s, "%.6lf", data->cur_time);
	strcat(s, "_plotina.dat");
	FILE* out = fopen(s, "w");

	fprintf(out, "VARIABLES = \"x\", \"ksi\", \"h\" \n");



		for (j = 0; j < n; ++j)
		{
			fprintf(out, "%lf %lf %lf\n", data->layer[1][j + 49 * n].x, data->layer[1][j + 49 * n].x / (data->cur_time*light_vel), data->layer[1][j + 49 * n].h);
		}
	
	fclose(out);
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void Print_layer(All_data *data, unsigned int l)
{
	char s[64];
	int i, j, k;
	double light_vel = sqrt(10. * 9.80665);
	unsigned int n = data->n, m = data->m;
	sprintf(s, "%.6lf", data->cur_time);
    strcat(s, ".dat");
	FILE* out = fopen(s, "w");

	fprintf(out, "FILETYPE = SOLUTION\n");
	fprintf(out, "VARIABLES = \"u\", \"v\", \"h\" \n");
	fprintf(out, "ZONE\nI=%d, J=%d\n ZONETYPE = Ordered, DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=CELLCENTERED)\n", n+1 , m+1 );
	fprintf(out, "SOLUTIONTIME=%lf\n", data->cur_time);
	
	for (i = 0; i < 3; ++i)
	{
		for (k = 0; k < m; ++k)
		{
			for (j = 0; j < n; ++j)
			{
				if ('y' == data->layer[1][i].square)
					{
						fprintf(out, "%lf ", 0.);
						continue;
					}
				if (0 == i)
					fprintf(out, "%lf ", data->layer[l][j + k*n].u);
				if (1 == i)
					fprintf(out, "%lf ", data->layer[l][j + k*n].v);
				if (2 == i)
					fprintf(out, "%lf ", data->layer[l][j + k*n].h);
			}
			fprintf(out, "\n");
		}
		fprintf(out, "\n");
	}
	fclose(out);
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void Init_data(All_data *data)
{
	int i;

	data->layer = (Point**) malloc(2 * sizeof(Point*));
	if (NULL == data->layer)
	{
		printf("data->layer malloc error\n");
		return ;
	}

	data->n = 200;
	data->m = 100;
	data->Fr = 1.;

	for (i = 0; i < 2; ++i)
	{
		data->layer[i] = (Point*) malloc(data->n*data->m*sizeof(Point));
		if (NULL == data->layer[i])
		{
			printf("data->layer[%d] malloc error\n", i);
			return ;
		}
	}
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void Init_net( All_data *data)
{
	int i;
	unsigned int j;
	unsigned int Jn, Jm;
	unsigned int n = data->n, m = data->m;
	data->d_x = 1. / n; data->d_y = data->d_x;
	double d_x = data->d_x, d_y = data->d_y;
	data->S = d_x*d_y;

	for (j = 0; j < m*n; ++j)
	{
		Jn = j%n;
		Jm = j/n;
		data->layer[0][j].x = Jn*d_x;
		data->layer[0][j].y = Jm*d_y;
		data->layer[1][j].x = data->layer[0][j].x;
		data->layer[1][j].y = data->layer[0][j].y;
		data->layer[0][j].square = 'n';
		data->layer[1][j].square = 'n';
		data->layer[0][j].c = 'n';
		data->layer[1][j].c = 'n';
	}
	Print_grid(data);
	//Create_Square(data, 140, 160);
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void Init_cond(All_data *data)
{
	int i;
	int Jm, Jn;
	unsigned int n = data->n, m = data->m;
	for (i = 0; i < m*n; ++i)
	{
		Jn = i % n;
		Jm = i / n;

		if ('y' == data->layer[1][i].square)
			continue;

		data->layer[0][i].h = 1.;
		data->layer[0][i].u = 0.;
		data->layer[0][i].v = 0.;

		data->layer[1][i].h = 1.;
		data->layer[1][i].u = 0.;
		data->layer[1][i].v = 0.;

		//walls
		if ((m - 1) == Jm)//up wall
		{
			data->layer[0][i].c = 'u';
			data->layer[1][i].c = 'u';
		}
		if (0 == Jm)//down wall
		{
			data->layer[0][i].c = 'd';
			data->layer[1][i].c = 'd';
		}
		if (0 == Jn)//left wall
		{
			data->layer[0][i].c = 'l';
			data->layer[1][i].c = 'l';
		}
		if ( (n-1) == Jn)//right wall
		{
			data->layer[0][i].c = 'r';
			data->layer[1][i].c = 'r';
		}
		
		if ( Jn < n/2 )
		{
			data->layer[0][i].h = 2.;
			data->layer[1][i].h = 2.;
		}
	}

	for (i = 0; i < m*n; ++i)
	{
		if ('y' == data->layer[1][i].square)
			continue;
		u2U(&data->layer[0][i]);
		u2U(&data->layer[1][i]);
	}
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void Stream_calculate(All_data *data, unsigned int j, struct Values1* val, char c )
{
	struct Values v1, v2, vb;
	double Fr = data->Fr;
	unsigned int n = data->n;
	switch (c)
	{
  	    case 'l':
	    {
					v2.rho = data->layer[0][j].h; v2.u = data->layer[0][j].u;
					v2.p = data->layer[0][j].h*data->layer[0][j].h / 2. / Fr / Fr;
					v1.rho = data->layer[0][j-1].h; v1.u = data->layer[0][j-1].u;
					v1.p = data->layer[0][j-1].h*data->layer[0][j-1].h / 2. / Fr / Fr;
					Crumbling(&v1, &v2, &vb, 0.);
					val->u = vb.u; 	val->h = vb.rho;
					val->v = data->layer[0][j].v;
					break;
	    }
	    case 'r':
	    {
					v1.rho = data->layer[0][j].h; v1.u = data->layer[0][j].u;
					v1.p = data->layer[0][j].h*data->layer[0][j].h / 2. / Fr / Fr;
					v2.rho = data->layer[0][j + 1].h; v2.u = data->layer[0][j + 1].u;
					v2.p = data->layer[0][j + 1].h*data->layer[0][j + 1].h / 2. / Fr / Fr;
					Crumbling(&v1, &v2, &vb, 0.);
					val->u = vb.u; 	val->h = vb.rho;
					val->v = data->layer[0][j].v;
					break;
	    }
	    case 'u':
	    {
					v1.rho = data->layer[0][j].h; v1.u = data->layer[0][j].v;
					v1.p = data->layer[0][j].h*data->layer[0][j].h / 2. / Fr / Fr;
					v2.rho = data->layer[0][j + n].h; v2.u = data->layer[0][j + n].v;
					v2.p = data->layer[0][j + n].h*data->layer[0][j + n].h / 2. / Fr / Fr;
					Crumbling(&v1, &v2, &vb, 0.);
					val->v = vb.u; 	val->h = vb.rho;
					val->u = data->layer[0][j].u;
					break;
	    }
	    case 'd':
	    {
					v2.rho = data->layer[0][j].h; v2.u = data->layer[0][j].v;
					v2.p = data->layer[0][j].h*data->layer[0][j].h / 2. / Fr / Fr;
					v1.rho = data->layer[0][j - n].h; v1.u = data->layer[0][j - n].v;
					v1.p = data->layer[0][j - n].h*data->layer[0][j - n].h / 2. / Fr / Fr;
					Crumbling(&v1, &v2, &vb, 0.);
					val->v = vb.u; 	val->h = vb.rho;
					val->u = data->layer[0][j].u;
					break;
	    }
	}


	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void SW1_Scheme(All_data *data, unsigned int j,
	double* Eu, double* Ed, double* Fl, double* Fr)
{
	double d_x = data->d_x;
	double d_t = data->d_t;
	double S = data->S;

	data->layer[1][j].u1 = data->layer[0][j].u1 - 
		(d_t*d_x / S)*(Eu[0] + Fr[0] - Ed[0] - Fl[0]);
	data->layer[1][j].u2 = data->layer[0][j].u2 - 
		(d_t*d_x / S)*(Eu[1] + Fr[1] - Ed[1] - Fl[1]);
	data->layer[1][j].u3 = data->layer[0][j].u3 - 
		(d_t*d_x / S)*(Eu[2] + Fr[2] - Ed[2] - Fl[2]);
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void Time_cycle(All_data *data)
{
	int i, counter = 0, point_check;
	int j;
	unsigned int Jn, Jm;
	unsigned int n = data->n, m = data->m;
	double max_vel = 0., temp_vel;
	double final_time = 3.;
	double Frud = data->Fr;
	double Fr2 = 1. / Frud / Frud;
	double Fl[3], Fr[3], Eu[3], Ed[3];
	Point* temp = NULL;
	struct Values1 *bu, *bd, *bl, *br;

	int start, now;
	start = clock();

	data->d_t = 0.00001;
	data->cur_time = 0.;

	while (data->cur_time < final_time)
	{
		/*refresh layer*/
		temp = data->layer[0];
		data->layer[0] = data->layer[1];
		data->layer[1] = temp;

		data->Q_l = (data->cur_time < 1.) ? data->cur_time : 1.;
		//data->Q_r = data->Q_l;
		//data->Q = 0.00001;

		/*Logic*/
		omp_set_num_threads(6);
#pragma omp parallel for private(Jn, Jm, bl, br, bu, bd, Eu, Ed, Fl, Fr)
		for (j = 0; j < n*m; ++j)
		{
			bu = (struct Values1 *) malloc(sizeof(struct Values1));
			bd = (struct Values1 *) malloc(sizeof(struct Values1));
			bl = (struct Values1 *) malloc(sizeof(struct Values1));
			br = (struct Values1 *) malloc(sizeof(struct Values1));
			if ('y' == data->layer[1][j].square)
				continue;

			Jn = j % n;
			Jm = j / n;
			{
				{
					//if (0 == Jn)/*left*/
					{
						//bl->h = data->layer[0][j].h;
						//bl->u = data->Q_l / bl->h;
						//bl->v = 0.;
					}
					//else 
					if (data->layer[0][j].c == 'l')/*left wall*/
					{
						bl->u = 0.;
						bl->v = data->layer[0][j].v;
						bl->h = data->layer[0][j].h;
					}
					else
						Stream_calculate(data, j, bl, 'l');/*calculating left stream*/
					B2F(bl, Fl, Frud);
				}

				{
				//if ((n - 1) == Jn)/*right*/
				{
					//br->h = 1.;
					//br->u = data->layer[1][j].h*data->layer[1][j].u;
					//br->v = data->layer[1][j].v;

				}
				//else
				if (data->layer[0][j].c == 'r')/*right wall*/
				{
					br->u = 0.;
					br->v = data->layer[0][j].v;
					br->h = data->layer[0][j].h;
				}
				else
					Stream_calculate(data, j, br, 'r');/*calculating right stream*/
				B2F(br, Fr, Frud);
			    }

				{
					if ((data->layer[0][j].c == 'u') || ((m - 1) == Jm))/*top wall*/
					{
						bu->u = data->layer[0][j].u;
						bu->v = 0.;
						bu->h = data->layer[0][j].h;
					}
					else
						Stream_calculate(data, j, bu, 'u');/*calculating top stream*/
					B2E(bu, Eu, Frud);
				}

				{
					if ((data->layer[0][j].c == 'd') || (0 == Jm))/*bottom wall*/
					{
						bd->u = data->layer[0][j].u;
						bd->v = 0.;
						bd->h = data->layer[0][j].h;
					}
					else
						Stream_calculate(data, j, bd, 'd');/*calculating down stream*/
					B2E(bd, Ed, Frud);
				}
			}
			SW1_Scheme(data, j, Eu, Ed, Fl, Fr);
			point_check = Point_check(&data->layer[1][j]);
			if (0 != point_check) /*NAN + h<0 check*/
			{
				printf("Problems in %d cell : %s\n", j, (point_check == 1)?"h<0":"NAN");
				free(bu); free(bd); free(bl); free(br);
				//return;
			}
			U2u(&data->layer[1][j]);
			
			free(bu); free(bd); free(bl); free(br);
		}

		/*Curant*/
		
		for (i = 0; i < n*m; ++i)
		{
			temp_vel = sqrt(data->layer[1][i].u*data->layer[1][i].u +
				data->layer[1][i].v*data->layer[1][i].v) + 
				sqrt(data->layer[1][i].h*Fr2);
			if (temp_vel > max_vel)
				max_vel = temp_vel;
		}
		
		data->d_t = (data->d_x / max_vel ) / 3.;
		data->cur_time += data->d_t;
		max_vel = 0.;

		/*Printing*/
		counter++;
		//Print_layer(data, 1);
		if (1==counter)
			Print_layer(data, 1);
		if (0 == (counter % 100))
		{
			now = clock();
			printf("%d iteration, current time = %lf, t=%lf\n", counter, data->cur_time, ((float)(now - start)) / CLOCKS_PER_SEC);
		}
		if (0 == (counter % 25))
		{
			//Print_layer(data, 1);
			Print_plotina(data);
		}
	}

	//free(bu); free(bd); free(bl); free(br);
	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
static void Del_data(All_data* data)
{
	int i;

	for (i = 0; i < 2; ++i)
		free(data->layer[i]);
	free(data->layer);

	return;
}
/*\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
int main(void)
{
	All_data data;

	Init_data(&data);
	Init_net(&data);
	Init_cond(&data);
	Time_cycle(&data);
	Del_data(&data);

	return 0;
}

