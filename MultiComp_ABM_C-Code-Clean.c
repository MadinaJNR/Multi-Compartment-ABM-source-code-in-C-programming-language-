/*================= Code changes and Notes =================================
This code is for the interaction of three agents MKK, ExR and MK. ExR are on the nuclear membrane and are in either an active or inactive states(with the ability to switch from one to another after certain number of iterations. The active state of ExR interacts with pMK (activated MK with state =2) while the inactive form doesn't do anything. After the interaction of pMK with the active ExR it translocates to th cytoplasm and in particular to a specific compartment in the cytoplasm. pMK will also change state from pMK to MK. MKK is found in compartment 1 and it's movemnet was restricted to that compartment. MKK interacted with MK and formed dMKK (55) and pMK which immediately transolcates to the nuclous.

Code address 1) translocation of MK into a special 3D compartement 
			 2) The movement restriction inside this 3D state of MKK and MK inside this 3D compartment
			 3) The interaction of MKK and MK (tests understanding of how the agents interact)
				A) The formation of both pMK and dMKK as a result of the interaction
				B) Transolcation of pMK to the nucleous as soon as it forms
			4) Interaction of pMK with active ExR (aExR) 
				A) To form MK which transolcates to its compartment of origion
				B) aExR deactivates and become dExR
			5) The activation of the 10 compartments 
			 
Changes made in comparision to MultiStateCode-ExR-MK+MK-MovementRestriction.c (aka MultiAgentCode-S5.3.c):

	ACTIVE compartments in the code [here is 10 comps]) and the (Y) referes to its being FINISHED and SUCCESSFUL:
	----> Addition of MKK and its three functions (inputdata, outputdata & move)
	----> Restricted its movement to one compartment to test the interaction with translocated MK and the formation of dMKK (55)	
	=====> The .xml file was changed to accomedate for MKK and addition of comptag hence the "i" in "Si" to accomedate for version changes in the .xml files >>>> FROM NOW ON CHANGES TO .XML WILL BE VERSIONED USING ROMAN NUMERALS 

Changes made in comparision to MultiAgentCode+MKK+2Compartments_v2.c:
	----> 1) Made MK from compa0 & compa1 transolcate to compa1
	----> 2) MK from compa1 & 0 move strictly in that compartment >>>> This IS NOT IMPLEMENTED IN MultiAgentCode-S5.1.c WHILE present IN S5.3.C
	----> 3) Changed the MKK_move function so that there is no compa0 for MKK >>> Compa0 is the nucleus!
	
Changes made in comparision to MultiAgentCode+MKK+2Compartments_v4.3V-3.c (In Dropbox named Y_MultiAgentCode-Cx10.c --> C referes to and x10 referes to the number:
	----> Just the activation of comp#3 and amking shure MKK & Mk behaved as wanted
	
----------------------------------------***General info***-----------------------------------------------
The .XML file accommedating code for this file has the same name with a .xml ending rather than .c.
	
The code for two compartments was finished by 12.04.2013 and was named MultiAgentCode+MKK+2Compartments_v4.3V-3.c 
	>>> Where (+) referes to an additional parameter or design of the code which is needed for its build up
	>>> Where _v(#) -> Referes to the versioin # of the code which achieved the purpose of the code 
	>>> Where v#.(#) referes to the sub-version of the working code which tackled the second issue of the code
	>>> Where v#.#-(Roman#) referes to subversion which tackled the 3rd issue of the code
	>>> Where v#.#-Roman#-(#) referes to subversion which improved an issue and tackeld an error

On 20.05.2013, the work was resumed on the code and the activation of the 3rd compartment was initiated >>> Hence the addition of Cx3
	>>> Cx3 refers to number of compartments (C) ---|> THIS WILL BE THE WAY TO USE IN THE NEW CODES FROM comp5
		>>> Cx3 takled the activation of MKK movment and its in comp#3 --|> Cx3(Roman#)
		>>> Cx3 tackled the activation of MK=>pMK ---|> Cx3ii-(#)
		>>> Cx3 tackled the restriction of MK movment in comp#3 ---|> Cx3ii-2(Roman#)

The name of the code addresses this as n in the naming ((Multi(...){n}compartments(...).c

FROM NOW ON CHANGES TO .XML WILL BE VERSIONED USING ROMAN NUMERALS 
	
	
=== Keys ===
=> = to
---> = These are (temp use for now [23.05.13])
===> = Because
=--> = Targeting
--|> = Resulting into (temp use for now [23.05.13])
# = number
>>> = Important to note the following statment in relation to the previous statment
-> = relates to

*/

#include "header.h"
#include <time.h>
#include <stdlib.h>

#define disrate_1 0 /* MT3 dissociation */

#define disrate_2 0 /* pMK dissociation to MKK and MK */

#define disrate_3 3/* MMT3 dissociation to MKK, trib and MK */						/* ### activated on 01062012 ### */

#define PI 3.141592
#define pi 3.141592
#define recdelay_constant 3
#define disK_constant 5

/* ExR code was finished on 26/11/2012 */

int MKK_outputdata()
{
	/* Send location message */
	add_MKKlocation_message(get_id(),		/* agent id */
						get_posx(),		/* agent x-axis position */
						get_posy(),		/* agent y-axis position */
						get_posz(),		/* agent z-axis position */
						get_state(),	/* agent type and state */
						get_iradius());	/* agent message range */
						
	
		
	/* if dormant MKK in the compartment was formed from MKK binding to MK */		
		if (get_state () == 55)
			{
			/* If it has a re-Activation Delay Period (RADP) value > 0 */
			 if( get_RADP() > 0)
				{
					// Then decrement RADP counter /
					set_RADP(get_RADP() - 1);
						// If processing has finished (RADP is zero) /
					if(get_RADP() == 0)
						{
						/* Change state from dormant MAPK to active MAPK and re-set the RADP value to a random value between zero and 100 [The value changes depending on the RADP value choosen and if the model used is stochastic or deterministic */
							set_state(0);
							set_RADP(rand()%(100-0)+0); 

						}
				
				}
			  }
			/* ====================================================================== */
					
	/* if dormant MKK [from interaction with Trib] in the compartment formed */
		if (get_state () == 555)
			{
			/* If it has a re-Activation Delay Period (RADP) value > 0 */
			 if( get_RADP() > 0)
				{
					// Then decrement RADP counter /
					set_RADP(get_RADP() - 1);
					// If processing has finished (RADP is zero) /
					if(get_RADP() == 0)
						{
							/* Change state from dormant MAPKK to active MAPKK and re-set the RADP value to 150 */
							set_RADP(150);
							set_state(0);
						}
				
				}
				
			}
		
		/* ====================================================================== */
			
	return 0;
}

int MKK_inputdata()
{
	xmachine_memory_MKK * xmemory = current_xmachine->xmachine_MKK;
	
	/* Create variables to use */
	double x1, y1, x2, y2, z1, z2, closestdist;
	int id2, s1, s2, closestid, closeststate; 
	double distance;
		
	/* Copy agent data to #1 variables */
	x1 = get_posx();
	y1 = get_posy();
	z1 = get_posz();
	s1 = get_state();
	
	/* Interaction radius of the agent (squared so don't need to sqrt compared variable) */
	closestdist = 99999 ;
	/* Index of closest agent to interact with and set to default as own id */
	closestid = get_id();

	/* Check for possible bindings */
	/* Get first location message */
	MKlocation_message = get_first_MKlocation_message();
	/* And loop through them until there are none left */
	

	while(MKlocation_message)
	{
		/* If the location message is not from the current agent */
		if((MKlocation_message->id != get_id()))
		{
			/* Copy location message data to #2 variables */
			id2 = MKlocation_message->id;
			x2 = MKlocation_message->x;
			y2 = MKlocation_message->y;
			z2 = MKlocation_message->z;
			s2 = MKlocation_message->state;
			
			/* Calculate the distance between agent and message sending agent */
			distance = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
			/* If distance is equal or smaller than the agent interaction radius */
			if(distance < closestdist)
			{
				
				/* if agent is free cytoplasmic MKK and the other is active MAPK */
				if((s1 == 0 && s2 == 200)) 				
				{
					/* if this agent is closer than the index of the closest */
					if(distance < closestdist)
					{
						/* update closest distance to bondable agent */
						closestdist = distance;
						/* also update the closest agent's id and state */
						closestid = id2;
						closeststate = s2;
					}
				}
			}
		}
		/* Get next location message */
		MKlocation_message = get_next_MKlocation_message(MKlocation_message);
	
	}
	
	/* If closest id is not the dummy value therefore a possible bond has been found */
	if(closestid != xmemory->id)
	{
		/* Send try to bond message */
		add_MKKnewbond_message(xmemory->id,         
						xmemory->state,
						closestid,           /* id of agent message is refering to */
						3,                   /* 3 is for a try to bond */
						closestdist,        /* distance to between agents */
						get_iradius(),
						get_posx(),
						get_posy(),
						get_posz());
	}	
	
	return 0;
}

int MKK_move()				
{
	xmachine_memory_MKK * xmemory = current_xmachine->xmachine_MKK;
	double cellradius = 10000.0;    	/* approximate cell radius in nano meters */
	double nuclradius = 4000.0;     	/* approximate cell nucleus radius in nano meters */
	double dt = 10.0;               	 /* Time differential */
	/* Molecule movement constants */
	double speed_ave = 2;				/* Average speed of protein agent */
	double speed_range = 1;		      /* Range of speed, + or - from average speed */
	double angle_range = pi/10.0;   /* Angle range in radians, + or - from current angle */
	/* Variables for calculations */
	double movex;
	double movey;
	double movez;

	/* Check if I have been bound */
	/* Get first bond message */
	MKfinalbond_message = get_first_MKfinalbond_message();
	/* And loop through messages until there are none left */
	while(MKfinalbond_message)
	{
		/* bond message refers to me (same id) and tag refers to a bind */
		if(MKfinalbond_message->idto == xmemory->id && MKfinalbond_message->bindunbind == 0)
		{
			/* Then check state and change accordingly */
			/* If MKK in cytoplasm and bound by MAPK; change state to dormant MKK [dMKK] */
			if(get_state() == 0  && MKfinalbond_message->statefrom == 200) 
				{
					set_state(55); 
				}
			/* If MKK in cytoplasm and bound by MAPK; change state to dormant MKK [dMKK] */
			if(xmemory->state == 0)
			{
				/* Remember ID of molecule bound to me */
				xmemory->boundindex = MKfinalbond_message->idfrom;
			}
		}
		/* Get next bond message */
		MKfinalbond_message = get_next_MKfinalbond_message(MKfinalbond_message);
	}
	/* If freely moving MAPKK molecule */
//	{
			if(xmemory->state==0 				/* MKK free in a compartment */
				|| xmemory->state==55			 /* dormant MKK (from MKK-MK) free in a compartment */
				|| xmemory->state==555)			/* dormant MKK (from MKK-Trib) free in a compartment */
			
			{
				/* If MKK agent reside within one of the cytoplasmic compartments */
				if(xmemory->comptag == 1
				|| xmemory->comptag == 2
				|| xmemory->comptag == 3
				|| xmemory->comptag == 4
				|| xmemory->comptag == 5
				|| xmemory->comptag == 6
				|| xmemory->comptag == 7
				|| xmemory->comptag == 8
				|| xmemory->comptag == 9
				|| xmemory->comptag == 10
				)
				{
			
						/* Calculate and update movement in polar coordinates */

						xmemory->movetheta = xmemory->movetheta + (-angle_range + 2.0*angle_range*((double)rand()/RAND_MAX));
						xmemory->movephi = xmemory->movephi + (-angle_range + 2.0*angle_range*((double)rand()/RAND_MAX));
						xmemory->mover = speed_ave + (-speed_range + 2.0*speed_range*((double)rand()/RAND_MAX));
						/* Calculate corresponding movement in cartesian coordinates */
						movex = xmemory->mover * cos(xmemory->movephi) * cos(xmemory->movetheta);
						movey = xmemory->mover * cos(xmemory->movephi) * sin(xmemory->movetheta);
						movez = xmemory->mover * sin(xmemory->movephi);
						/* Update cartesian positions with respect to time */
						xmemory->posx += movex * dt;
						xmemory->posy += movey * dt;
						xmemory->posz += movez * dt;
						/* Calculate polar position from cartesian position */
						xmemory->postheta = atan2(xmemory->posy, xmemory->posx);
						xmemory->posphi = atan2(xmemory->posz,sqrt((xmemory->posx*xmemory->posx)+(xmemory->posy*xmemory->posy)));
						xmemory->posr = sqrt((xmemory->posx*xmemory->posx)+(xmemory->posy*xmemory->posy)+(xmemory->posz*xmemory->posz));	
					
												
					/* =========== MKK movement restriction to within the compartment MKK reside in =================== */

			/* ------------------------------------------------------------------------------------------------------------------------------------ */		
						/* If MKK reside within compartment 1 */
						if(xmemory->comptag == 1)
						{
								/* If  7000 < MKK's X-coordinate < 8000 and past compartment 1; update X-coordinate */
								if( xmemory->posx > 8000 || xmemory->posx < 7000) 
								
								 {
								 /* If MKK x-coordinate > 8000 */
									if( xmemory->posx > 8000)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = 8000 -(xmemory->posx-8000); 
									}
									/* If MKK x-coordinate < 7000 */
									if( xmemory->posx < 7000)
									{
									/* mirror MKK into the compartment */
										xmemory->posx =  (7000 - xmemory->posx ) + 7000; 
									}
								  }	
								/* If  4248 < MKK's Y-coordinate < 5248 and past compartment 1; update Y-coordinate */	
								else if( xmemory->posy > 5248 || xmemory->posy < 4248)
								
								{
									 /* When MKK y-coordinate > 5248 */
									if( xmemory->posy > 5248)
									{
									/* mirror MKK into the compartment */
										xmemory->posy = (5248 - xmemory->posy) + 5248; 
									}
									 /* When MKK x-coordinate < 4248 */
									if( xmemory->posy < 4248)
									{
									/* mirror MKK into the compartment */	
										xmemory->posy = (4248 - xmemory->posy) + 4248; 
									}
								 }

								/* If  -7200 > the MKK's Z-coordinate > -6000 and past compartment 1; update Z-coordinate */
								else if( xmemory->posz < -7200 || xmemory->posz > -6000)
								{
									 /* When MKK z-coordinate < -7200 */								
									if( xmemory->posz < -7200)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (-7200 - xmemory->posz) - 7200; 
									}
									/* When MKK z-coordinate > -6000 */								
									if( xmemory->posz > -6000)
									{
										/* mirror MKK into the compartment */
										xmemory->posz =  (-6000 - xmemory->posz) - 6000; 
									}
								}
						}	
				
			/* ------------------------------------------------------------------------------------------------------------------------------------ */		
							/* If MKK reside within compartment 2 */
							if(xmemory->comptag == 2)	
							{
										/* If  6800 < MKK's X-coordinate < 5800 and past compartment 2; update X-coordinate */				
										if( xmemory->posx > 6800 || xmemory->posx < 5800) 
										{
											/* When MKK x-coordinate > 6800 */										
											if( xmemory->posx > 6800)
											{
												/* mirror MKK into the compartment */
												xmemory->posx = 6800 -(xmemory->posx-6800); 
											}
											/* When MKK x-coordinate < 5800 */										
											if( xmemory->posx < 5800)
											{
												/* mirror MKK into the compartment */
												xmemory->posx =  (5800 - xmemory->posx ) + 5800; 
											}
										}
										/* If  6048 < MKK's Y-coordinate < 5048 and past compartment 2; update Y-coordinate */		
										else if( xmemory->posy > 6048 || xmemory->posy < 5048)
										{
											/* When MKK y-coordinate > 6048 */									
											if( xmemory->posy > 6048)
											{
												xmemory->posy = (6048 - xmemory->posy) + 6048; 
											}
											/* When MKK y-coordinate < 5048 */	
											if( xmemory->posy < 5048)
											{
												/* mirror MKK into the compartment */
												xmemory->posy = (5048 - xmemory->posy) + 5048; 
											}
										 }
										/* If  7923 > the MKK's Z-coordinate < 6923 and past compartment 2; update Z-coordinate */
										else if( xmemory->posz < 7923 || xmemory->posz > 6923)
										{
											/* When MKK z-coordinate > 79232*/												
											if( xmemory->posz > 7923)
											{
												/* mirror MKK into the compartment */
												xmemory->posz = (7923 - xmemory->posz) + 7923; 
											}
											/* When MKK z-coordinate < 6923 */		
											if( xmemory->posz < 6923 )
											{
												/* mirror MKK into the compartment */
												xmemory->posz =  (6923 - xmemory->posz)+ 6923; 
											}
										}
							}
			/* ------------------------------------------------------------------------------------------------------------------------------------ */		

							/* If MKK reside within compartment 3 */
							if(xmemory->comptag == 3)
							{
								/* If  1700 < MKK's X-coordinate < 900 and past compartment 3; update X-coordinate */				
								if( xmemory->posx > 1700 || xmemory->posx < 900) 
								{
									/* When MKK x-coordinate > 1700 */
									if( xmemory->posx > 1700)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = 1700 -(xmemory->posx-1700); 
									}
									/* When MKK x-coordinate < 900 */
									if( xmemory->posx < 900)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (900 - xmemory->posx ) + 900; 
									}
								 }	
								/* If  4200 < MKK's Y-coordinate < 3600 and past compartment 3; update Y-coordinate */			
								else if( xmemory->posy > 4200 || xmemory->posy < 3600)
								{
									/* When MKK y-coordinate > 4200 */	
									if( xmemory->posy > 4200)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4200 - xmemory->posy) + 4200; 
									}
									/* When MKK y-coordinate < 3600 */	
									if( xmemory->posy < 3600)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (3600 - xmemory->posy) + 3600; 
									}
								 }

								/* If  3600 > the MKK's Z-coordinate < 4600 and past compartment 3; update Z-coordinate */
								else if( xmemory->posz < 3600 || xmemory->posz > 4350)
								{
									/* When MKK z-coordinate > 4600*/	
									if( xmemory->posz > 4350)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (4350 - xmemory->posz) + 4350; 
									}
									/* When MKK z-coordinate < 3600 */	
									if( xmemory->posz < 3600)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (3600 - xmemory->posz) + 3600; 
									}
								}
							}			
							
			/* ------------------------------------------------------------------------------------------------------------------------------------ */		
							/* If MKK reside within compartment 4 */
							if(xmemory->comptag == 4)
							{
								/* If  4400 < MKK's X-coordinate < 3700 and past compartment 4; update X-coordinate */				
								if( xmemory->posx > 4400 || xmemory->posx < 3700) 
								{
									/* When MKK x-coordinate > 4400 */
									if( xmemory->posx > 4400)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = 4400 -(xmemory->posx-4400); 
									}
									/* When MKK x-coordinate < 3700 */
									if( xmemory->posx < 3700)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (3700 - xmemory->posx ) + 3700; 
									}
								 }	
									
								/* If   5580 < MKK's Y-coordinate < 4900 and past compartment 4; update Y-coordinate */			
								else if( xmemory->posy > 5580 || xmemory->posy < 4900)
								{
									/* When MKK y-coordinate > 5580 */	
									if( xmemory->posy > 5580)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (5580 - xmemory->posy) + 5580; 
									}
									
									/* When MKK y-coordinate < 4900 */	
									if( xmemory->posy < 4900)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4900 - xmemory->posy) + 4900; 
									}
								 }

								
								/* If  900 > the MKK's Z-coordinate < 100 and past compartment 4; update Z-coordinate */
								else if( xmemory->posz < 100 || xmemory->posz > 900)
								{
									/* When MKK z-coordinate < 900 */
									if( xmemory->posz > 900)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (900 - xmemory->posz) + 900; 
									}
									
									/* When MKK z-coordinate < 100 */
									if( xmemory->posz < 100)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (100 - xmemory->posz) + 100; 
									}
								}
							}
			/* ------------------------------------------------------------------------------------------------------------------------------------ */
							/* If MKK reside within compartment 5 */
							if(xmemory->comptag == 5)
							{
								/* If  -2550 < MKK's X-coordinate < 3380 and past compartment 5; update X-coordinate */				
							  if( xmemory->posx > -2550 || xmemory->posx < -3380) 
								
								 {
									/* When MKK x-coordinate > -2550 */
									if( xmemory->posx > -2550)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = -2550 -(xmemory->posx - -2550 ); 
									}
									/* When MKK x-coordinate < -3380 */
									if( xmemory->posx < -3380)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (-3380 - xmemory->posx )  + -3380; 
									}
								  }
								 /* If   -2550 < MKK's Y-coordinate < -1730 and past compartment 5; update Y-coordinate */			
								else if( xmemory->posy > -1730 || xmemory->posy < -2550) 
								{
									/* When MKK y-coordinate > -1730 */	
									if( xmemory->posy > -1730)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = -1730 -(xmemory->posy - -1730 ); 
									}
									
									/* When MKK y-coordinate < -2500 */	 
									if( xmemory->posy < -2550)
									{
										/* mirror MKK into the compartment */
										xmemory->posy =  (-2550 - xmemory->posy )  + -2550; 
									}
								  }
								/* If  -5130 > the MKK's Z-coordinate < -5950 and past compartment 5; update Z-coordinate */
								else if( xmemory->posz > -5130 || xmemory->posz < -5950) 
								{
									/* When MKK z-coordinate > -5130 */
									if( xmemory->posz > -5130)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = -5130 -(xmemory->posz + 5130 ); 
									}
									/* When MKK z-coordinate < -5950*/
									if( xmemory->posz < -5950)
									{
										/* mirror MKK into the compartment */
										xmemory->posz =  -5950 - (xmemory->posz  + 5950); 
									}
								}
								
							}
			/* ------------------------------------------------------------------------------------------------------------------------------------ */
								/* If MKK reside within compartment 6 */
								if(xmemory->comptag == 6)
								{
									/* If  -5380 < MKK's X-coordinate < -4550 and past compartment 6; update X-coordinate */				
									if( xmemory->posx > -4550 || xmemory->posx < -5380) 
									
									 {
										/* When MKK x-coordinate > -4550 */										
										if( xmemory->posx > -4550)
										{
											/* mirror MKK into the compartment */
											xmemory->posx = -4550 -(xmemory->posx - -4550 ); 
										}
										/* When MKK x-coordinate > -5380 */
										if( xmemory->posx < -5380)
										{
											/* mirror MKK into the compartment */
											xmemory->posx =  (-5380 - xmemory->posx )  + -5380; 
										}
									  }	
										
									/* If   4580 < MKK's Y-coordinate < 5400 and past compartment 6; update Y-coordinate */			
									else if( xmemory->posy > 5400 || xmemory->posy < 4580)
									{
										/* When MKK y-coordinate > 5400 */	 
										if( xmemory->posy > 5400)
										{
											/* mirror MKK into the compartment */
											xmemory->posy = (5400 - xmemory->posy)  + 5400; 
										}
										
										/* When MKK y-coordinate < -4580 */	 
										if( xmemory->posy < 4580)
										{
											/* mirror MKK into the compartment */
											xmemory->posy = (4580 - xmemory->posy)  + 4580; 
										}
									 }
 
									
									/* If  -5130 > the MKK's Z-coordinate < -5950 and past compartment 6; update Z-coordinate */
									else if( xmemory->posz > -5130 || xmemory->posz < -5950) 
									{
										/* When MKK z-coordinate > -5130*/
										if( xmemory->posz > -5130)
										{
											/* mirror MKK into the compartment */
											xmemory->posz = -5130 -(xmemory->posz + 5130 ); 
										}
										
										/* When MKK z-coordinate < -5950*/
										if( xmemory->posz < -5950)
										{
											/* mirror MKK into the compartment */
											xmemory->posz =  -5950 - (xmemory->posz  + 5950); 
										}
									}	
																
							}
							
			/* ------------------------------------------------------------------------------------------------------------------------------------ */
								/* If MKK reside within compartment 7 */
								if(xmemory->comptag == 7)
								{
									/* If  -3650 < MKK's X-coordinate < -4480 and past compartment 7; update X-coordinate */				
									if( xmemory->posx > -3650 || xmemory->posx < -4480) 
									
									 {
										/* When MKK x-coordinate > -3650 */
										if( xmemory->posx > -3650)
										{
											/* mirror MKK into the compartment */
											xmemory->posx = -3650 -(xmemory->posx - -3650 ); 
										}
										/* When MKK x-coordinate < -4480 */
										if( xmemory->posx < -4480)
										{
											/* mirror MKK into the compartment */ 
											xmemory->posx =  (-4480 - xmemory->posx )  + -4480; 
										}
									  }	
										
									/* If   -5500 < MKK's Y-coordinate < -6300 and past compartment 7; update Y-coordinate */			
									else if( xmemory->posy > -5500 || xmemory->posy < -6300) 
									{
										/* When MKK y-coordinate > -5500*/	 
										if( xmemory->posy > -5500)
										{
											/* mirror MKK into the compartment */
											xmemory->posy = -5500 -(xmemory->posy - -5500 ); 
										}
										/* When MKK y-coordinate < -6300*/	 
										if( xmemory->posy < -6300)
										{
											/* mirror MKK into the compartment */
											xmemory->posy =  (-6300 - xmemory->posy )  + -6300; 
										}
																										
									}			
									
								/* If  6300 > the MKK's Z-coordinate < 5450 and past compartment 7; update Z-coordinate */
								else if( xmemory->posz > 6300 || xmemory->posz < 5450)
								{
									/* When MKK z-coordinate < 6300 */
									if( xmemory->posz > 6300)
									{
										/* mirror MKK into the compartment */
										xmemory->posz =  5887; 
									}
									
									/* When MKK z-coordinate < 5450*/
									if( xmemory->posz < 5450)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = 5890; 
									}
								
								}
									
					
							 }
							
							
			/* ------------------------------------------------------------------------------------------------------------------------------------ */
							/* If MKK reside within compartment 8 */
							if(xmemory->comptag == 8)
							{
								/* If  5600 < MKK's X-coordinate < 6250 and past compartment 8; update X-coordinate */				
								if( xmemory->posx > 6250 || xmemory->posx < 5600) 
								
								 {
									/* When MKK x-coordinate > 6250 */
									if( xmemory->posx > 6250)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = 6250 -(xmemory->posx - 6250 ); 
									}
									/* When MKK x-coordinate <5600 */
									if( xmemory->posx < 5600)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (5600 - xmemory->posx )  + 5600; 
									}
								  }	
								/* If   4200 < MKK's Y-coordinate < 4880 and past compartment 8; update Y-coordinate */				
								else if( xmemory->posy > 4880 || xmemory->posy < 4200)
								{
									/* When MKK y-coordinate > 4880 */	 
									if( xmemory->posy > 4880)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4880 - xmemory->posy)  + 4880; 
									}
									/* When MKK y-coordinate < 4200 */	 
									if( xmemory->posy < 4200)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4200 - xmemory->posy)  + 4200; 
									}
								}

								/* If  4700 < the MKK's Z-coordinate < 3900 and past compartment 8; update Z-coordinate */
								else if( xmemory->posz < 3900 || xmemory->posz > 4700)
								{
									/* When MKK z-coordinate > 4700 */
									if( xmemory->posz > 4700)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (4700 - xmemory->posz)  + 4700; 
									}
									/* When MKK z-coordinate < 3900*/
									if( xmemory->posz < 3900)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (3900 - xmemory->posz)  + 3900 ; 
									}
								}
							}
			
			
			
			/* ------------------------------------------------------------------------------------------------------------------------------------ */		
					/* If MKK reside within compartment 9 */
					if(xmemory->comptag == 9)
							{
								/* If  -4650 < MKK's X-coordinate < -5350 and past compartment 9; update X-coordinate */		
								if( xmemory->posx > -4650 || xmemory->posx < -5350) 
								
								 {
									/* When MKK x-coordinate > -4650 */
									if( xmemory->posx > -4650)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = -4650 -(xmemory->posx - -4650 ); 
									}
									/* When MKK x-coordinate < -5350 */
									if( xmemory->posx < -5350)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (-5350 - xmemory->posx )  + -5350; 
									}
								  }	
									
								/* If   5560 < MKK's Y-coordinate < 4850 and past compartment 9; update Y-coordinate */				
								else if( xmemory->posy > 5560 || xmemory->posy < 4850)
								{
									/* When MKK y-coordinate > 5560 */
									if( xmemory->posy > 5560)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (5560 - xmemory->posy)  + 5560; 
									}
									/* When MKK y-coordinate < 4850 */
									if( xmemory->posy < 4850)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4850 - xmemory->posy)  + 4860; 
									}
								 }

								/* If  6300 < the MKK's Z-coordinate < 5590 and past compartment 9; update Z-coordinate */
								else if( xmemory->posz < 5590 || xmemory->posz > 6300)
								{
									/* When MKK z-coordinate > 6300*/
									if( xmemory->posz > 6300)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (6300 - xmemory->posz)  + 6300; 
									}
									/* When MKK z-coordinate < 5590*/
									if( xmemory->posz < 5590)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (5590 - xmemory->posz)  + 5590 ; 
									}
								}
							}
			
			
			/* ------------------------------------------------------------------------------------------------------------------------------------ */		
							/* If MKK reside within compartment 10 */
							if(xmemory->comptag == 10)
							{
								/* If  5880 < MKK's X-coordinate < 5170 and past compartment 10; update X-coordinate */
								if( xmemory->posx > 5880 || xmemory->posx < 5170) 
								 {
									/* When MKK x-coordinate < 5880 */
									if( xmemory->posx > 5880)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = (5880 - xmemory->posx) + 5880 ; 
									}
									/* When MKK x-coordinate < 5170 */
									if( xmemory->posx < 5170)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (5170 - xmemory->posx )  + 5170; 
									}
								  }	
									
								/* If   5530 < MKK's Y-coordinate < 4750 and past compartment 10; update Y-coordinate */				
								else if( xmemory->posy > 5530 || xmemory->posy < 4750)
								{
									/* When MKK y-coordinate > 5530 */
									if( xmemory->posy > 5530)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (5530 - xmemory->posy)  + 5530; 
									}
									/* When MKK y-coordinate < 4750 */
									if( xmemory->posy < 4750)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4750 - xmemory->posy)  + 4750; 
									}
								 }

								/* If  -5480 < the MKK's Z-coordinate < -4800 and past compartment 10; update Z-coordinate */
								else if( xmemory->posz > -4800 || xmemory->posz < -5480) 
								 {
									/* When MKK z-coordinate < -4800*/
									if( xmemory->posz > -4800)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = -4800 -(xmemory->posz + 4800 ); 
									}
									/* When MKK z-coordinate < -5480*/
									if( xmemory->posz < -5480)
									{
										/* mirror MKK into the compartment */
										xmemory->posz =  -5480 - (xmemory->posz  + 5480); 
									}
								  }	
																
							}
				/* ------------------------------------------------------------------------------------------------------------------------------------ */		
							
										
					}
				
				}
						

			
//	}
	return 0;
}

int MK_outputdata()
{
	/* Send location message */
	add_MKlocation_message(get_id(),		/* agent id */
						get_posx(),		/* agent x-axis position */
						get_posy(),		/* agent y-axis position */
						get_posz(),		/* agent z-axis position */
						get_state(),	/* agent type and state */
						get_iradius());	/* agent message range */
											
		/* ====================================================================== */
		
	/* If protein agent state is dormant gene-expression activator (DGEA) phospho-MAPK */			
		if (get_state () == 24)
			{
			/* If it has a RADP8 value > 0 */
			 if( get_RADP8() > 0)
				{
					// Then decrement RADP8 counter /
					set_RADP8(get_RADP8() - 1);
					// If processing has finished (RADP is zero) /
					if(get_RADP8() == 0)
						{
							// Change state from DGEA pMAPK to active pMAPK */
							set_state(2);
						}
				
				}
			}

			/* ====================================================================== */
					
	/* if MK state is dormant phospho-MAPK (dpMK) in the nucleus */		
		if (get_state () == 42)
			{
			 if( get_RADP8() > 0)
				{
					// Then decrement RADP8 counter /
					set_RADP8(get_RADP8() - 1);
					// If processing has finished (RADP is zero) /
					if(get_RADP8() == 0)
						{
							// Change state from DGEA pMAPK to active pMAPK */
							set_state(2);
						}
				
				}
				
			}
		
		/* ====================================================================== */
			
	return 0;
}

int MK_inputdata()
{
	/* Create variables to use */
	double x1, y1, x2, y2, z1, z2, closestdist;
	int id2, s1, s2, closestid, closeststate; 
	double distance;
		
	/* Copy agent data to #1 variables */
	x1 = get_posx();
	y1 = get_posy();
	z1 = get_posz();
	s1 = get_state();
	
	/* Interaction radius of the agent (squared so don't need to sqrt compared variable) */
	closestdist = 99999 ;//get_iradius()*get_iradius();
	/* Index of closest agent to interact with and set to default as own id */
	closestid = get_id();

	
		/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

	/* Check for possible bindings */
	/* Get first location message */
	ExRlocation_message = get_first_ExRlocation_message();
	/* And loop through them until there are none left */
	

	while(ExRlocation_message)
	{
		/* If the location message is not from the current agent */
		if((ExRlocation_message->id != get_id()))
		{
			/* Copy location message data to #2 variables */
			id2 = ExRlocation_message->id;
			x2 = ExRlocation_message->x;
			y2 = ExRlocation_message->y;
			z2 = ExRlocation_message->z;
			s2 = ExRlocation_message->state;
			
			/* Calculate the distance between agent and message sending agent */
			/* Do not square root as more efficient to square interaction radius */
			distance = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
			/* If distance is equal or smaller than the agent interaction radius */
			if(distance < closestdist)
			{
				if(
				/* if agent is free nuclear pMK and the other is active and/or dormant Exrec*/
				(s1 == 2 && s2 == 501) 				
				|| (s1 == 2 && s2 == 500)  

				)

				{
		 
					/* if this agent is closer than the index of the closest */
					if(distance < closestdist)
					{
						/* update closest distance to bondable agent */
						closestdist = distance;
						/* also update the closest agent's id and state */
						closestid = id2;
						closeststate = s2;
					}
				}
			}
		}
		/* Get next location message */
		ExRlocation_message = get_next_ExRlocation_message(ExRlocation_message);
		
	}
	
	while(MKKlocation_message) 
	{
		/* If the location message is not from the current agent */
		if((MKKlocation_message->id != get_id()))
		{
			/* Copy location message data to #2 variables */
			id2 = MKKlocation_message->id;
			x2 = MKKlocation_message->x;
			y2 = MKKlocation_message->y;
			z2 = MKKlocation_message->z;
			s2 = MKKlocation_message->state;
			
			/* Calculate the distance between agent and message sending agent */
			/* Do not square root as more efficient to square interaction radius */
			distance = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
			/* If distance is equal or smaller than the agent interaction radius */
			if(distance < closestdist)
			{
				if(
				/* if agent is free cytoplasmic MAPK and the other active MKK */
				(s1 == 200 && s2 == 0) 				
				)

				{
					/* if this agent is closer than the index of the closest */
					if(distance < closestdist)
					{
						/* update closest distance to bondable agent */
						closestdist = distance;
						/* also update the closest agent's id and state */
						closestid = id2;
						closeststate = s2;
					}
				}
			}
		}
		/* Get next location message */
		MKKlocation_message = get_next_MKKlocation_message(MKKlocation_message);
	}	
	
	return 0;
}
int MK_checkbondtries()
{
	double cellradius = 10000.0;    /* approximate cell radius in nano meters */
	double nuclradius = 4000.0;     /* approximate cell nucleus radius in nano meters */
	
	xmachine_memory_MK * xmemory = current_xmachine->xmachine_MK;
	
	double closestdist = 0.0;
	int closestid = xmemory->id; // dummy value
	int closeststate = 0;

	/* Check if I have been bound to ExR */
	/* Get first bond message from ExR agent*/
	ExRnewbond_message = get_first_ExRnewbond_message();
	/* And loop through messages until there are none left */
	while(ExRnewbond_message)
	{
		if(ExRnewbond_message->idto == get_id() && ExRnewbond_message->bindunbind == 3)
		{
			if((closestid == xmemory->id) || (ExRnewbond_message->distance < closestdist))
			{
				closestdist = ExRnewbond_message->distance;
				closestid = ExRnewbond_message->idfrom;
				closeststate = ExRnewbond_message->statefrom;
			}
		}
		
		/* Get next trybond message */
		ExRnewbond_message = get_next_ExRnewbond_message(ExRnewbond_message);
	}
	
	/* Check if I have been bound to MKK */
	/* Get first bond message from MKK agent*/
	MKKnewbond_message = get_first_MKKnewbond_message();
	/* And loop through messages until there are none left */
	while(MKKnewbond_message)
	{
		if(MKKnewbond_message->idto == get_id() && MKKnewbond_message->bindunbind == 3)
		{
			if((closestid == xmemory->id) || (MKKnewbond_message->distance < closestdist))
			{
				closestdist = MKKnewbond_message->distance;
				closestid = MKKnewbond_message->idfrom;
				closeststate = MKKnewbond_message->statefrom;
			}
		}
		
		/* Get next trybond message */
		MKKnewbond_message = get_next_MKKnewbond_message(MKKnewbond_message);
	}
	
	// if not dummy value then bond created
	if(closestid != xmemory->id)
	{
		/* Send bond message to molecule to bind */
		add_MKfinalbond_message(xmemory->id,
						xmemory->state,
						closestid,          /* id of MKK agent message is refering to */
						0.0,    /* distance not used */
						get_iradius(),
						get_posx(),
						get_posy(),
						get_posz());
		
		// ***************************************************************  (START #) ******************************************************************************** \\ 
						/* If pMK free in nucleus and bound by exporting nuclear receptor then make pMK bound to exporting nuclear receptor [ExR], change state to active MAPK, re-set RADP and followed by export out of the nucleus */
						if(get_state() == 2  && closeststate == 500) 
						{
							set_state(200);											
							set_RADP3(90);
	/* =============================== MK translocation into the corresponding compartment ============================================================= */	
		
		/* 		-------------------------------------------------------------------------------------------------------------------------------------------------------		*/			
								/* if MK originated from compartment 10, translocate back into the compartment with randomised 3D coordinates  */ 
								if (xmemory->comptag == 10)
								{	
									set_posx(rand()%(5880-5170)+5170); 
									set_posy(rand()%(5530-4750)+4750); 
									set_posz(rand()%(-4780+5550)-5550);
								}
								
								/* if MK originated from compartment 9, translocate back into the compartment #9 with randomised 3D coordinates  */ 
								if (xmemory->comptag == 9)
								{	
									set_posx(rand()%(-4650+5350)-5350); 
									set_posy(rand()%(5560-4850)+4850); 
									set_posz(rand()%(6300-5590)+5590);	
								}
		
								/* if MK originated from compartment 8, translocate back into the compartment #8 with randomised 3D coordinates  */ 
								if (xmemory->comptag == 8)
								{	
									set_posx(rand()%(6250-5600)+5600); 
									set_posy(rand()%(4880-4200)+4200); 
									set_posz(rand()%(4750-4000)+4000);	
								}
		
								/* if MK originated from compartment 7, translocate back into the compartment #7 with randomised 3D coordinates  */ 
								if (xmemory->comptag == 7)
								{	
									set_posx(rand()%(-4480+3650)-3650);
									set_posy(rand()%(-6300+5500)-5500);
									set_posz(rand()%(6300-5450)+5450); 
								}
		
								/* if MK originated from compartment 6, translocate back into the compartment #6 with randomised 3D coordinates  */ 
								if (xmemory->comptag == 6)
								{	
									set_posx(rand()%(-4550+5380)-5380);
									set_posy(rand()%(5400-4580)+4580); 
									set_posz(rand()%(-5130+5950)-5950);	
								}
		
								/* if MK originated from compartment 5, translocate back into the compartment #5 with randomised 3D coordinates  */ 
								if (xmemory->comptag == 5)
								{	
									set_posx(rand()%(-2550+3380)-3380); 
									set_posz(rand()%(-1730+2550)-2550);	
									set_posz(rand()%(-5130+5950)-5950);	
								}
		
								/* if MK originated from compartment 4, translocate back into the compartment #4 with randomised 3D coordinates  */ 
								if (xmemory->comptag == 4)
								{	
									set_posx(rand()%(4400-3700)+3700); 
									set_posy(rand()%(5580-4900)+4900); 
									set_posz(rand()%(900-100)+100);	
								}
								
								/* if MK originated from compartment 3, translocate back into the compartment #3 with randomised 3D coordinates  */ 
								if (xmemory->comptag == 3)
								{	
									set_posx(rand()%(1700-900)+900); 
									set_posy(rand()%(4200-3600)+3600); 
									set_posz(rand()%(4380-3600)+3600);	
								}
	
								/* if MK originated from compartment 2, translocate back into the compartment #2 with randomised 3D coordinates  */ 
								if (xmemory->comptag == 2)
								{	
									set_posx(rand()%(6800-5800)+5800); 
									set_posy(rand()%(6048-5048)+5048); 
									set_posz(rand()%(7923-6923)+6923);	
								}	
							
							/* if MK originated from compartment 1 and the nucleus [compartment 0], translocate back into the compartment #1 with randomised 3D coordinates  */ 
								if(xmemory->comptag == 0 || xmemory->comptag == 1)
								{
									set_posx(rand()%(8000-7000)+7000);
									set_posy(rand()%(5248-4248)+4248);
									
									set_posz(rand()%(8123-7123)+-7123);		
								
								
								}
					

							}	
	/* =============================== MK translocation into the corresponding compartment ============================================================= */	
	
						/* If MAPK free in the cytoplasm and bound to MAPKK change state to pMK and transolcate into the nucleous regardless of which compartment MK originated from*/	
						if(get_state() == 200  && closeststate == 0) 
						{
							if(xmemory->comptag == 2
							|| xmemory->comptag == 1
							|| xmemory->comptag == 0 
							|| xmemory->comptag == 3
							|| xmemory->comptag == 4
							|| xmemory->comptag == 5
							|| xmemory->comptag == 6
							|| xmemory->comptag == 7
							|| xmemory->comptag == 8
							|| xmemory->comptag == 9
							|| xmemory->comptag == 10
							) 
							{
							set_state(2);
							set_mover(0.0);
							}
							 
					
						}
	}

	return 0;
}

int MK_move()				/* **** where i got on 27 & 28.11.2012 **** */
{
	xmachine_memory_MK * xmemory = current_xmachine->xmachine_MK;
	double cellradius = 10000.0;    /* approximate cell radius in nano meters */
	double nuclradius = 4000.0;     /* approximate cell nucleus radius in nano meters */
	double dt = 10.0;                /* Time differential */
	/* Molecule movement constants */
	double speed_ave = 2;//0.0;        /* Average speed of protein agent */
	double speed_range = 1;//0.0;      /* Range of speed, + or - from average speed */
	double angle_range = pi/10.0;   /* Angle range in radians, + or - from current angle */
	/* Variables for calculations */
	double movex;
	double movey;
	double movez;


	/* If freely moving pMAPK */
	if(( xmemory->state==2          /* pMK free in nucleous*/																
		|| xmemory->state==24			 /* dormant pMK free in nucleus */
		|| xmemory->state==42			 /* dormant pMK free in nucleus */
		))
				
	{
		/* Calculate and update movement in polar coordinates */

		xmemory->movetheta = xmemory->movetheta + (-angle_range + 2.0*angle_range*((double)rand()/RAND_MAX));
		xmemory->movephi = xmemory->movephi + (-angle_range + 2.0*angle_range*((double)rand()/RAND_MAX));
		xmemory->mover = speed_ave + (-speed_range + 2.0*speed_range*((double)rand()/RAND_MAX));
		/* Calculate corresponding movement in cartesian coordinates */
		movex = xmemory->mover * cos(xmemory->movephi) * cos(xmemory->movetheta);
		movey = xmemory->mover * cos(xmemory->movephi) * sin(xmemory->movetheta);
		movez = xmemory->mover * sin(xmemory->movephi);
		/* Update cartesian positions with respect to time */
		xmemory->posx += movex * dt;
		xmemory->posy += movey * dt;
		xmemory->posz += movez * dt;
		/* Calculate polar position from cartesian position */
		xmemory->postheta = atan2(xmemory->posy, xmemory->posx);
		xmemory->posphi = atan2(xmemory->posz,sqrt((xmemory->posx*xmemory->posx)+(xmemory->posy*xmemory->posy)));
		xmemory->posr = sqrt((xmemory->posx*xmemory->posx)+(xmemory->posy*xmemory->posy)+(xmemory->posz*xmemory->posz));
		
		/* If position is beyond the nuclear membrane */
				
				if(xmemory->posr > nuclradius)
				{
					/* Mirror position back into nucleus */
					xmemory->posr = (2*nuclradius)-xmemory->posr;
					/* Update new position */
					xmemory->posx = xmemory->posr * cos(xmemory->posphi) * cos(xmemory->postheta);
					xmemory->posy = xmemory->posr * cos(xmemory->posphi) * sin(xmemory->postheta);
					xmemory->posz = xmemory->posr * sin(xmemory->posphi);
					/* Mirror direction of movement */
					movex = -movex;
					movey = -movey;
					movez = -movez;
					/* Update new movement */
					xmemory->movetheta = atan2(movey, movex);
					xmemory->movephi = atan2(movez,sqrt((movex*movex)+(movey*movey)));
					xmemory->mover = sqrt((movex*movex)+(movey*movey)+(movez*movez));
				}
				
	}
			/* If MAPK and resides whithin a cytoplasmic compartment */
			if(xmemory->state==200) /* If in cytoplasm */

			{
					if(xmemory->comptag == 1
					|| xmemory->comptag == 0 
					|| xmemory->comptag == 2
					|| xmemory->comptag == 3
					|| xmemory->comptag == 4
					|| xmemory->comptag == 5
					|| xmemory->comptag == 6
					|| xmemory->comptag == 7
					|| xmemory->comptag == 8
					|| xmemory->comptag == 9
					|| xmemory->comptag == 10
					)
					
					{	/* Calculate and update movement in polar coordinates */

						xmemory->movetheta = xmemory->movetheta + (-angle_range + 2.0*angle_range*((double)rand()/RAND_MAX));
						xmemory->movephi = xmemory->movephi + (-angle_range + 2.0*angle_range*((double)rand()/RAND_MAX));
						xmemory->mover = speed_ave + (-speed_range + 2.0*speed_range*((double)rand()/RAND_MAX));
						/* Calculate corresponding movement in cartesian coordinates */
						movex = xmemory->mover * cos(xmemory->movephi) * cos(xmemory->movetheta);
						movey = xmemory->mover * cos(xmemory->movephi) * sin(xmemory->movetheta);
						movez = xmemory->mover * sin(xmemory->movephi);
						/* Update cartesian positions with respect to time */
						xmemory->posx += movex * dt;
						xmemory->posy += movey * dt;
						xmemory->posz += movez * dt;
						/* Calculate polar position from cartesian position */
						xmemory->postheta = atan2(xmemory->posy, xmemory->posx);
						xmemory->posphi = atan2(xmemory->posz,sqrt((xmemory->posx*xmemory->posx)+(xmemory->posy*xmemory->posy)));
						xmemory->posr = sqrt((xmemory->posx*xmemory->posx)+(xmemory->posy*xmemory->posy)+(xmemory->posz*xmemory->posz));	
						
					/* Calculate redirection of molecules if leaving prescribed compartment */
					
					
				/* Update new position */

							/* =========== MKK movement restriction to within the compartment MKK reside in =================== */

					/* ------------------------------------------------------------------------------------------------------------------------------------ */		

						if(xmemory->comptag == 0 || xmemory->comptag == 1)
						{
								/* If  7000 < MK's X-coordinate < 8000 and past compartment 1; update X-coordinate */
								if( xmemory->posx > 8000 || xmemory->posx < 7000) 
								
								 {
								 /* If MK x-coordinate > 8000 */
									if( xmemory->posx > 8000)
									{
										/* mirror MK into the compartment */
										xmemory->posx = 8000 -(xmemory->posx-8000); 
									}
									/* If MK x-coordinate < 7000 */
									if( xmemory->posx < 7000)
									{
									/* mirror MK into the compartment */
										xmemory->posx =  (7000 - xmemory->posx ) + 7000; 
									}
								  }	
								/* If  4248 < MK's Y-coordinate < 5248 and past compartment 1; update Y-coordinate */	
								else if( xmemory->posy > 5248 || xmemory->posy < 4248)
								
								{
									 /* When MK y-coordinate > 5248 */
									if( xmemory->posy > 5248)
									{
									/* mirror MK into the compartment */
										xmemory->posy = (5248 - xmemory->posy) + 5248; 
									}
									 /* When MK x-coordinate < 4248 */
									if( xmemory->posy < 4248)
									{
									/* mirror MK into the compartment */	
										xmemory->posy = (4248 - xmemory->posy) + 4248; 
									}
								 }

								/* If  -7200 > the MK's Z-coordinate > -6000 and past compartment 1; update Z-coordinate */
								else if( xmemory->posz < -7200 || xmemory->posz > -6000)
								{
									 /* When MK z-coordinate < -7200 */								
									if( xmemory->posz < -7200)
									{
										/* mirror MK into the compartment */
										xmemory->posz = (-7200 - xmemory->posz) - 7200; 
									}
									/* When MK z-coordinate > -6000 */								
									if( xmemory->posz > -6000)
									{
										/* mirror MK into the compartment */
										xmemory->posz =  (-6000 - xmemory->posz) - 6000; 
									}
								}
					 }

					/* ------------------------------------------------------------------------------------------------------------------------------------ */		

						/* If MKK reside within compartment 2 */
							if(xmemory->comptag == 2)	
							{
										/* If  6800 < MKK's X-coordinate < 5800 and past compartment 2; update X-coordinate */				
										if( xmemory->posx > 6800 || xmemory->posx < 5800) 
										{
											/* When MKK x-coordinate > 6800 */										
											if( xmemory->posx > 6800)
											{
												/* mirror MKK into the compartment */
												xmemory->posx = 6800 -(xmemory->posx-6800); 
											}
											/* When MKK x-coordinate < 5800 */										
											if( xmemory->posx < 5800)
											{
												/* mirror MKK into the compartment */
												xmemory->posx =  (5800 - xmemory->posx ) + 5800; 
											}
										}
										/* If  6048 < MKK's Y-coordinate < 5048 and past compartment 2; update Y-coordinate */		
										else if( xmemory->posy > 6048 || xmemory->posy < 5048)
										{
											/* When MKK y-coordinate > 6048 */									
											if( xmemory->posy > 6048)
											{
												xmemory->posy = (6048 - xmemory->posy) + 6048; 
											}
											/* When MKK y-coordinate < 5048 */	
											if( xmemory->posy < 5048)
											{
												/* mirror MKK into the compartment */
												xmemory->posy = (5048 - xmemory->posy) + 5048; 
											}
										 }
										/* If  7923 > the MKK's Z-coordinate < 6923 and past compartment 2; update Z-coordinate */
										else if( xmemory->posz < 7923 || xmemory->posz > 6923)
										{
											/* When MKK z-coordinate > 79232*/												
											if( xmemory->posz > 7923)
											{
												/* mirror MKK into the compartment */
												xmemory->posz = (7923 - xmemory->posz) + 7923; 
											}
											/* When MKK z-coordinate < 6923 */		
											if( xmemory->posz < 6923 )
											{
												/* mirror MKK into the compartment */
												xmemory->posz =  (6923 - xmemory->posz)+ 6923; 
											}
										}
							}
					/* ------------------------------------------------------------------------------------------------------------------------------------ */		

			
							if(xmemory->comptag == 3)
							{
								/* If  1700 < MKK's X-coordinate < 900 and past compartment 3; update X-coordinate */				
								if( xmemory->posx > 1700 || xmemory->posx < 900) 
								{
									/* When MKK x-coordinate > 1700 */
									if( xmemory->posx > 1700)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = 1700 -(xmemory->posx-1700); 
									}
									/* When MKK x-coordinate < 900 */
									if( xmemory->posx < 900)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (900 - xmemory->posx ) + 900; 
									}
								 }	
								/* If  4200 < MKK's Y-coordinate < 3600 and past compartment 3; update Y-coordinate */			
								else if( xmemory->posy > 4200 || xmemory->posy < 3600)
								{
									/* When MKK y-coordinate > 4200 */	
									if( xmemory->posy > 4200)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4200 - xmemory->posy) + 4200; 
									}
									/* When MKK y-coordinate < 3600 */	
									if( xmemory->posy < 3600)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (3600 - xmemory->posy) + 3600; 
									}
								 }

								/* If  3600 > the MKK's Z-coordinate < 4600 and past compartment 3; update Z-coordinate */
								else if( xmemory->posz < 3600 || xmemory->posz > 4350)
								{
									/* When MKK z-coordinate > 4600*/	
									if( xmemory->posz > 4350)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (4350 - xmemory->posz) + 4350; 
									}
									/* When MKK z-coordinate < 3600 */	
									if( xmemory->posz < 3600)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (3600 - xmemory->posz) + 3600; 
									}
								}
							}			
							
			/* ------------------------------------------------------------------------------------------------------------------------------------ */
							/* If MKK reside within compartment 4 */
							if(xmemory->comptag == 4)
							{
								/* If  4400 < MKK's X-coordinate < 3700 and past compartment 4; update X-coordinate */				
								if( xmemory->posx > 4400 || xmemory->posx < 3700) 
								{
									/* When MKK x-coordinate > 4400 */
									if( xmemory->posx > 4400)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = 4400 -(xmemory->posx-4400); 
									}
									/* When MKK x-coordinate < 3700 */
									if( xmemory->posx < 3700)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (3700 - xmemory->posx ) + 3700; 
									}
								 }	
									
								/* If   5580 < MKK's Y-coordinate < 4900 and past compartment 4; update Y-coordinate */			
								else if( xmemory->posy > 5580 || xmemory->posy < 4900)
								{
									/* When MKK y-coordinate > 5580 */	
									if( xmemory->posy > 5580)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (5580 - xmemory->posy) + 5580; 
									}
									
									/* When MKK y-coordinate < 4900 */	
									if( xmemory->posy < 4900)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4900 - xmemory->posy) + 4900; 
									}
								 }

								
								/* If  900 > the MKK's Z-coordinate < 100 and past compartment 4; update Z-coordinate */
								else if( xmemory->posz < 100 || xmemory->posz > 900)
								{
									/* When MKK z-coordinate < 900 */
									if( xmemory->posz > 900)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (900 - xmemory->posz) + 900; 
									}
									
									/* When MKK z-coordinate < 100 */
									if( xmemory->posz < 100)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (100 - xmemory->posz) + 100; 
									}
								}
							}
			/* ------------------------------------------------------------------------------------------------------------------------------------ */
							
							/* If MKK reside within compartment 5 */
							if(xmemory->comptag == 5)
							{
								/* If  -2550 < MKK's X-coordinate < 3380 and past compartment 5; update X-coordinate */				
							  if( xmemory->posx > -2550 || xmemory->posx < -3380) 
								
								 {
									/* When MKK x-coordinate > -2550 */
									if( xmemory->posx > -2550)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = -2550 -(xmemory->posx - -2550 ); 
									}
									/* When MKK x-coordinate < -3380 */
									if( xmemory->posx < -3380)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (-3380 - xmemory->posx )  + -3380; 
									}
								  }
								 /* If   -2550 < MKK's Y-coordinate < -1730 and past compartment 5; update Y-coordinate */			
								else if( xmemory->posy > -1730 || xmemory->posy < -2550) 
								{
									/* When MKK y-coordinate > -1730 */	
									if( xmemory->posy > -1730)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = -1730 -(xmemory->posy - -1730 ); 
									}
									
									/* When MKK y-coordinate < -2500 */	 
									if( xmemory->posy < -2550)
									{
										/* mirror MKK into the compartment */
										xmemory->posy =  (-2550 - xmemory->posy )  + -2550; 
									}
								  }
								/* If  -5130 > the MKK's Z-coordinate < -5950 and past compartment 5; update Z-coordinate */
								else if( xmemory->posz > -5130 || xmemory->posz < -5950) 
								{
									/* When MKK z-coordinate > -5130 */
									if( xmemory->posz > -5130)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = -5130 -(xmemory->posz + 5130 ); 
									}
									/* When MKK z-coordinate < -5950*/
									if( xmemory->posz < -5950)
									{
										/* mirror MKK into the compartment */
										xmemory->posz =  -5950 - (xmemory->posz  + 5950); 
									}
								}
								
							}
							
			/* ------------------------------------------------------------------------------------------------------------------------------------ */
								/* If MKK reside within compartment 6 */
								if(xmemory->comptag == 6)
								{
									/* If  -5380 < MKK's X-coordinate < -4550 and past compartment 6; update X-coordinate */				
									if( xmemory->posx > -4550 || xmemory->posx < -5380) 
									
									 {
										/* When MKK x-coordinate > -4550 */										
										if( xmemory->posx > -4550)
										{
											/* mirror MKK into the compartment */
											xmemory->posx = -4550 -(xmemory->posx - -4550 ); 
										}
										/* When MKK x-coordinate > -5380 */
										if( xmemory->posx < -5380)
										{
											/* mirror MKK into the compartment */
											xmemory->posx =  (-5380 - xmemory->posx )  + -5380; 
										}
									  }	
										
									/* If   4580 < MKK's Y-coordinate < 5400 and past compartment 6; update Y-coordinate */			
									else if( xmemory->posy > 5400 || xmemory->posy < 4580)
									{
										/* When MKK y-coordinate > 5400 */	 
										if( xmemory->posy > 5400)
										{
											/* mirror MKK into the compartment */
											xmemory->posy = (5400 - xmemory->posy)  + 5400; 
										}
										
										/* When MKK y-coordinate < -4580 */	 
										if( xmemory->posy < 4580)
										{
											/* mirror MKK into the compartment */
											xmemory->posy = (4580 - xmemory->posy)  + 4580; 
										}
									 }
 
									
									/* If  -5130 > the MKK's Z-coordinate < -5950 and past compartment 6; update Z-coordinate */
									else if( xmemory->posz > -5130 || xmemory->posz < -5950) 
									{
										/* When MKK z-coordinate > -5130*/
										if( xmemory->posz > -5130)
										{
											/* mirror MKK into the compartment */
											xmemory->posz = -5130 -(xmemory->posz + 5130 ); 
										}
										
										/* When MKK z-coordinate < -5950*/
										if( xmemory->posz < -5950)
										{
											/* mirror MKK into the compartment */
											xmemory->posz =  -5950 - (xmemory->posz  + 5950); 
										}
									}	
																
							}
							
				/* ------------------------------------------------------------------------------------------------------------------------------------ */
								
								/* If MKK reside within compartment 7 */
								if(xmemory->comptag == 7)
								{
									/* If  -3650 < MKK's X-coordinate < -4480 and past compartment 7; update X-coordinate */				
									if( xmemory->posx > -3650 || xmemory->posx < -4480) 
									
									 {
										/* When MKK x-coordinate > -3650 */
										if( xmemory->posx > -3650)
										{
											/* mirror MKK into the compartment */
											xmemory->posx = -3650 -(xmemory->posx - -3650 ); 
										}
										/* When MKK x-coordinate < -4480 */
										if( xmemory->posx < -4480)
										{
											/* mirror MKK into the compartment */ 
											xmemory->posx =  (-4480 - xmemory->posx )  + -4480; 
										}
									  }	
										
									/* If   -5500 < MKK's Y-coordinate < -6300 and past compartment 7; update Y-coordinate */			
									else if( xmemory->posy > -5500 || xmemory->posy < -6300) 
									{
										/* When MKK y-coordinate > -5500*/	 
										if( xmemory->posy > -5500)
										{
											/* mirror MKK into the compartment */
											xmemory->posy = -5500 -(xmemory->posy - -5500 ); 
										}
										/* When MKK y-coordinate < -6300*/	 
										if( xmemory->posy < -6300)
										{
											/* mirror MKK into the compartment */
											xmemory->posy =  (-6300 - xmemory->posy )  + -6300; 
										}
																										
									}			
									
									/* If  6300 > the MKK's Z-coordinate < 5450 and past compartment 7; update Z-coordinate */
									else if( xmemory->posz > 6300 || xmemory->posz < 5450)
									{
										/* When MKK z-coordinate < 6300 */
										if( xmemory->posz > 6300)
										{
											/* mirror MKK into the compartment */
											xmemory->posz =  5887; 
										}
										
										/* When MKK z-coordinate < 5450*/
										if( xmemory->posz < 5450)
										{
											/* mirror MKK into the compartment */
											xmemory->posz = 5890; 
										}
									
									}
										
						
								 }
							
			/* ------------------------------------------------------------------------------------------------------------------------------------ */
			
							/* If MKK reside within compartment 8 */
							if(xmemory->comptag == 8)
							{
								/* If  5600 < MKK's X-coordinate < 6250 and past compartment 8; update X-coordinate */				
								if( xmemory->posx > 6250 || xmemory->posx < 5600) 
								
								 {
									/* When MKK x-coordinate > 6250 */
									if( xmemory->posx > 6250)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = 6250 -(xmemory->posx - 6250 ); 
									}
									/* When MKK x-coordinate <5600 */
									if( xmemory->posx < 5600)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (5600 - xmemory->posx )  + 5600; 
									}
								  }	
								/* If   4200 < MKK's Y-coordinate < 4880 and past compartment 8; update Y-coordinate */				
								else if( xmemory->posy > 4880 || xmemory->posy < 4200)
								{
									/* When MKK y-coordinate > 4880 */	 
									if( xmemory->posy > 4880)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4880 - xmemory->posy)  + 4880; 
									}
									/* When MKK y-coordinate < 4200 */	 
									if( xmemory->posy < 4200)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4200 - xmemory->posy)  + 4200; 
									}
								}

								/* If  4700 < the MKK's Z-coordinate < 3900 and past compartment 8; update Z-coordinate */
								else if( xmemory->posz < 3900 || xmemory->posz > 4700)
								{
									/* When MKK z-coordinate > 4700 */
									if( xmemory->posz > 4700)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (4700 - xmemory->posz)  + 4700; 
									}
									/* When MKK z-coordinate < 3900*/
									if( xmemory->posz < 3900)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (3900 - xmemory->posz)  + 3900 ; 
									}
								}
							}
			
			
			
			/* ------------------------------------------------------------------------------------------------------------------------------------ */		
						/* If MKK reside within compartment 9 */
							if(xmemory->comptag == 9)
							{
								/* If  -4650 < MKK's X-coordinate < -5350 and past compartment 9; update X-coordinate */		
								if( xmemory->posx > -4650 || xmemory->posx < -5350) 
								
								 {
									/* When MKK x-coordinate > -4650 */
									if( xmemory->posx > -4650)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = -4650 -(xmemory->posx - -4650 ); 
									}
									/* When MKK x-coordinate < -5350 */
									if( xmemory->posx < -5350)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (-5350 - xmemory->posx )  + -5350; 
									}
								  }	
									
								/* If   5560 < MKK's Y-coordinate < 4850 and past compartment 9; update Y-coordinate */				
								else if( xmemory->posy > 5560 || xmemory->posy < 4850)
								{
									/* When MKK y-coordinate > 5560 */
									if( xmemory->posy > 5560)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (5560 - xmemory->posy)  + 5560; 
									}
									/* When MKK y-coordinate < 4850 */
									if( xmemory->posy < 4850)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4850 - xmemory->posy)  + 4860; 
									}
								 }

								/* If  6300 < the MKK's Z-coordinate < 5590 and past compartment 9; update Z-coordinate */
								else if( xmemory->posz < 5590 || xmemory->posz > 6300)
								{
									/* When MKK z-coordinate > 6300*/
									if( xmemory->posz > 6300)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (6300 - xmemory->posz)  + 6300; 
									}
									/* When MKK z-coordinate < 5590*/
									if( xmemory->posz < 5590)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = (5590 - xmemory->posz)  + 5590 ; 
									}
								}
							}
/* ------------------------------------------------------------------------------------------------------------------------------------ */			
			
							/* If MKK reside within compartment 10 */
							if(xmemory->comptag == 10)
							{
								/* If  5880 < MKK's X-coordinate < 5170 and past compartment 10; update X-coordinate */
								if( xmemory->posx > 5880 || xmemory->posx < 5170) 
								 {
									/* When MKK x-coordinate < 5880 */
									if( xmemory->posx > 5880)
									{
										/* mirror MKK into the compartment */
										xmemory->posx = (5880 - xmemory->posx) + 5880 ; 
									}
									/* When MKK x-coordinate < 5170 */
									if( xmemory->posx < 5170)
									{
										/* mirror MKK into the compartment */
										xmemory->posx =  (5170 - xmemory->posx )  + 5170; 
									}
								  }	
									
								/* If   5530 < MKK's Y-coordinate < 4750 and past compartment 10; update Y-coordinate */				
								else if( xmemory->posy > 5530 || xmemory->posy < 4750)
								{
									/* When MKK y-coordinate > 5530 */
									if( xmemory->posy > 5530)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (5530 - xmemory->posy)  + 5530; 
									}
									/* When MKK y-coordinate < 4750 */
									if( xmemory->posy < 4750)
									{
										/* mirror MKK into the compartment */
										xmemory->posy = (4750 - xmemory->posy)  + 4750; 
									}
								 }

								/* If  -5480 < the MKK's Z-coordinate < -4800 and past compartment 10; update Z-coordinate */
								else if( xmemory->posz > -4800 || xmemory->posz < -5480) 
								 {
									/* When MKK z-coordinate < -4800*/
									if( xmemory->posz > -4800)
									{
										/* mirror MKK into the compartment */
										xmemory->posz = -4800 -(xmemory->posz + 4800 ); 
									}
									/* When MKK z-coordinate < -5480*/
									if( xmemory->posz < -5480)
									{
										/* mirror MKK into the compartment */
										xmemory->posz =  -5480 - (xmemory->posz  + 5480); 
									}
								  }	
																
							}
			/* ------------------------------------------------------------------------------------------------------------------------------------ */								  
			
					}


					
				}
	return 0;
}

int ExR_outputdata() 
{
	/* Send location message */
	add_ExRlocation_message(get_id(),		/* agent id */
						get_posx(),		/* agent x-axis position */
						get_posy(),		/* agent y-axis position */
						get_posz(),		/* agent z-axis position */
						get_state(),	/* agent type and state */
						get_iradius());	/* agent message range */


		if(get_recdelay() > 0)
		{
			/* Then decrement receptor delay counter */
			set_recdelay(get_recdelay() - 1);
			/* If processing has finished (receptor delay is zero) */
			if(get_recdelay() == 0)
			{
				/* If inactive ExR [dExR], revert back to active ExR */
				if(get_state() == 501) set_state(500);
				/* If active ExR, revert back to inactive ExR [dExR] */
				else if(get_state() == 500) set_state(501);
				
			}
		}
		/* If recdelay is zero, re-set redelay to 10*/
		if(get_recdelay() == 0) set_recdelay(10);

	
	return 0;
}

int ExR_inputdata()
{
	xmachine_memory_ExR * xmemory = current_xmachine->xmachine_ExR;
	
	/* Create variables to use */
	double x1, y1, x2, y2, z1, z2, closestdist;
	int id2, s1, s2, closestid, closeststate;
	double distance;
	
	/* Copy agent data to #1 variables */
	x1 = get_posx();
	y1 = get_posy();
	z1 = get_posz();
	s1 = get_state();
	/* Interaction radius of the agent (squared so don't need to sqrt compared var) */
	closestdist = get_iradius()*get_iradius();
	/* Index of closest agent to interact with and set to default as own id */
	closestid = get_id();

	/* Check for possible bindings */
	/* Get first location message */
	MKlocation_message = get_first_MKlocation_message();
	/* And loop through them until there are none left */
	while(MKlocation_message)
	{
		/* If the location message is not from the current agent */
		if((MKlocation_message->id != get_id()))
		{
			/* Copy location message data to #2 variables */
			id2 = MKlocation_message->id;
			x2 = MKlocation_message->x;
			y2 = MKlocation_message->y;
			z2 = MKlocation_message->z;
			s2 = MKlocation_message->state;
			
			/* Calculate the distance between agent and message sending agent */
			/* Do not square root as more efficient to square interaction radius */
			distance = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
			/* If distance is equal or smaller than the agent interaction radius */
			if(distance < closestdist)
			{
				/* If agent is active Exporting Nuclear Receptor [ExR] */
				if((s1 == 500 && s2 == 2))
					
				{
					/* if this agent is closer than the index of the closest */
					if(distance < closestdist)
					{
						/* update closest distance to bondable agent */
						closestdist = distance;
						/* also update the closest agent's id and state */
						closestid = id2;
						closeststate = s2;
					}
				}
			}
			
		}
		/* Get next location message */
		MKlocation_message = get_next_MKlocation_message(MKlocation_message);
	}
		
	/* If closest id is not the dummy value therefore a possible bond has been found */
	if(closestid != xmemory->id)
	{
		/* Send try to bond message */
		add_ExRnewbond_message(xmemory->id,         
						xmemory->state,
						closestid,           /* id of agent message is refering to */
						3,                   /* 3 is for a try to bond */
						closestdist,        /* distance to between agents */
						get_iradius(),
						get_posx(),
						get_posy(),
						get_posz());
	}
	
	return 0;
}
int ExR_move()
{
	xmachine_memory_ExR * xmemory = current_xmachine->xmachine_ExR;
	
	double nuclradius = 4000.0;     /* Approximate cell nucleus radius in nano meters */		
	double dt = 10.0;                /* Time differential */									
	/* Nuclear receptor constants */
	double nuc_speed_ave = 2.0;               /* Average speed of  ... */ 
	double nuc_speed_range = 1.0;                /* Range of speed, + or - from average speed */
	double nuc_ang_dev = pi/10.0;               	/* Angle dev */  
	double nuc_ang_speed_ave = nuc_speed_ave;    	/* Angle speed average */
	double nuc_ang_speed_dev = nuc_speed_range;  	/* Angle speed dev */
	
	/* Check if I have been bound */
	/* Get first bond message */
	MKfinalbond_message = get_first_MKfinalbond_message();
	/* And loop through messages until there are none left */
	while(MKfinalbond_message)
	{
		/* bond message refers to me (same id) and tag refers to a bind */
		if(MKfinalbond_message->idto == xmemory->id && MKfinalbond_message->bindunbind == 0)
		{

		/* If agent is active Exporting Nuclear Receptor [ExR] bound to pMAPK; change state to dormant ExR and re-set recdelay to 20*/
	if(get_state() == 500  && MKfinalbond_message->statefrom == 2) 
						{
						set_state(501); 
						set_recdelay (20);
						}
			/* If active exporting receptor */
			if(xmemory->state == 500)
			{
				/* Remember ID of molecule bound to me */
				xmemory->boundindex = MKfinalbond_message->idfrom;
			}
			/* If dorment exporting receptor */
			if(xmemory->state == 501)
			{
				/* Remember ID of molecule bound to me */
				xmemory->boundindex = MKfinalbond_message->idfrom;
			}
		}
		/* Get next bond message */
		MKfinalbond_message = get_next_MKfinalbond_message(MKfinalbond_message);
	}

	/* If nuclear receptor */
	if(xmemory->state >= 500 && xmemory->state < 610)
	{
		/* Calculate movement */
		double ang_change = -nuc_ang_dev + 2.0*((double)rand()/RAND_MAX);
		double theta_temp = xmemory->postheta + xmemory->movetheta + ang_change;
		double phi_temp = xmemory->posphi + xmemory->movephi - ang_change;
		double mov_nucrec = nuc_ang_speed_ave - nuc_ang_speed_dev + 2.0*nuc_ang_speed_dev*((double)rand()/RAND_MAX);
		/* Update position variables in memory */
		xmemory->postheta = xmemory->postheta + dt*mov_nucrec*(theta_temp - xmemory->postheta);
		xmemory->posphi = xmemory->posphi + dt*mov_nucrec*(phi_temp - xmemory->posphi);
		xmemory->posx = xmemory->posr * cos(xmemory->posphi) * cos(xmemory->postheta);
		xmemory->posy = xmemory->posr * cos(xmemory->posphi) * sin(xmemory->postheta);
		xmemory->posz = xmemory->posr * sin(xmemory->posphi);
	}
	
	return 0;
}
