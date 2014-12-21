#include "header.h"
#include <time.h>
#include <stdlib.h>

#define disrate_1 0 /* MT3 dissociation */

#define disrate_2 0 /* pMK dissociation to MKK and MK */

#define disrate_3 3/* MMT3 dissociation to MKK, trib and MK */						/* ### activated on 01062012 ### */

#define PI 3.141592
#define pi 3.141592
#define recdelay_constant 3
#define RADP_constant 5

int Protein_outputdata()
{
	/* Send location message */
	add_location_message(get_id(),		/* Protein agent id */
						get_posx(),		/* Protein agent x-axis position */
						get_posy(),		/* Protein agent y-axis position */
						get_posz(),		/* Protein agent z-axis position */
						get_state(),	/* Protein agent type and state */
						get_iradius());	/* Protein agent message range */
	
																
	/* If protein agent state is dormant MAPK which interacted with Trib protein agent */	
	if (get_state () == 222)
			{
			/* If it has a re-Activation Delay Period3 (RADP3) value > 0 */
			 if( get_RADP3() > 0 )
				{
					// Then decrement RADP counter /
					set_RADP3(get_RADP3() - 1);
					// If processing has finished (RADP is zero) /
					if(get_RADP3() == 0)
						{
							/* Change state from dormant MAPK to active MAPK */
							set_state(2);
						}
				
				}
			}
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
							// Set RADP to 9 and change state from DGEA pMAPK to active pMAPK */
							set_RADP8(9);
							set_state(2);
						}
				
				}
			}
		/* If protein agent state is dormant phospho-MAPK */			
		if (get_state () == 42)
			{
			 if( get_RADP8() > 0)	
				{
					// Then decrement RADP counter /
					set_RADP8(get_RADP8() - 1);
					// If processing has finished (receptor delay is zero) /
					if(get_RADP8() == 0)
						{
							// Send bond message to bound molecule to unbind /
							set_RADP8(8);
							set_state(2);
						}
				
				}
			}
			
		if (get_state () == 21)
			{
			 if( get_RADP5() > 0)	
				{
					// Then decrement RADP counter /
					set_RADP5(get_RADP5() - 1);
					// If processing has finished (receptor delay is zero) /
					if(get_RADP5() == 0)
						{
							// Send bond message to bound molecule to unbind /
							set_RADP5(5);
							set_state(400);
						}
				
				}
			}
			
		if (get_state () == 22)
			{
			 if( get_RADP5() > 0)	
				{
					// Then decrement RADP counter /
					set_RADP5(get_RADP5() - 1);
					// If processing has finished (receptor delay is zero) /
					if(get_RADP5() == 0)
						{
							// Send bond message to bound molecule to unbind /
							set_RADP5(5);
							set_state(403);
						}
				
				}
			}
	/* If protein agent state is Trib which interacted with protein agent state MAPKK*/			
	if (get_state () == 1)
			{
			/* If it has a re-Activation Delay Period2 (RADP2) value > 0 */
			 if( get_RADP2() > 0)
				{
					// Then decrement RADP counter /
					set_RADP2(get_RADP2() - 1);
					// If processing has finished (RADP is zero) /
					if(get_RADP2() == 0)
						{
							/* Change state from reacted-Trib to active Trib */
							set_state(900);
						}
				
				}
			}		
		/* If protein agent state is Trib which interacted with protein agent state MAPK */			
	if (get_state () == 11)
			{
			/* If it has a re-Activation Delay Period4 (RADP4) value > 0*/
			 if( get_RADP4() > 0)
				{
					// Then decrement RADP4 counter /
					set_RADP4(get_RADP4() - 1);
					// If processing has finished (RADP4 is zero) /
					if(get_RADP4() == 0)
						{
							// Set RADP2 to 25 and revert state from MAPK reacted-Trib to MAPKK reacted-Trib //
							set_RADP2(25);
							set_state(1);
						}
				
				}
			}	
		
	return 0;
}

int Protein_inputdata()
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
	/* Interaction radius of the agent */
	closestdist = 99999 ;
	/* Index of closest agent to interact with and set to default as own id */
	closestid = get_id();
	
	/* If it has a re-Activation Delay Period (RADP) value > 0*/					
	if(get_RADP() > 0)
		{
			/* Then decrement RADP counter */
			set_RADP(get_RADP() - 1);
			/* If processing has finished (RADP is zero) */
			if(get_RADP() == 0)
			{
				/* Revert dMAPKK back to MAPKK */
				if(get_state() == 55) set_state(0);
			}
		}
			/* If processing has finished (RADP is zero) re-set RADP to 10 */
		if(get_RADP() == 0) set_RADP(10);
	
	/* If it has a re-Activation Delay Period6 (RADP6) value > 0*/	
	if(get_RADP6() > 0)
		{
			/* Then decrement RADP6 counter */
			set_RADP6(get_RADP6() - 1);
			/* If processing has finished6 (RADP6 is zero) */
			if(get_RADP6() == 0)
			{
				/* Revert SB-TF back to SF-TF */
				 if(get_state() == 400) set_state(401); 
				 /* Revert SF-TF back to SB-TF */
				else if(get_state() == 401) set_state(400);
			}
		}
		/* If processing has finished (RADP6 is zero) re-set RADP6 to 5 */
		if(get_RADP6() == 0) set_RADP6(5);
	
/* If it has a re-Activation Delay Period6 (RADP6) value > 0*/	
	if(get_RADP7() > 0)
		{
			/* Then decrement RADP7 counter */
			set_RADP7(get_RADP7() - 1);
			/* If processing has finished7 (RADP7 is zero) */
			if(get_RADP7() == 0)
			{
				/* Revert MB-TF back to MF-TF */
				 if(get_state() == 403) set_state(404);
				 /* Revert MF-TF back to MB-TF */
				else if(get_state() == 404) set_state(403);
			}
		}
		/* If processing has finished (RADP6 is zero) re-set RADP6 to 5 */
		if(get_RADP7() == 0) set_RADP7(5);
	
	 /* Check if I have been unbound */
	/* Get first bond message */
	startbond_message = get_first_startbond_message();
	/* And loop through them until there are none left */
	while(startbond_message)
	{
		/* Get next bond message */
		startbond_message = get_next_startbond_message(startbond_message);
	}

	/* Check for possible bindings */
	/* Get first location message */
	location_message = get_first_location_message();
	/* And loop through them until there are none left */
	

	while(location_message)
	{
		/* If the location message is not from the current agent */
		if((location_message->id != get_id()))
		{
			/* Copy location message data to #2 variables */
			id2 = location_message->id;
			x2 = location_message->x;
			y2 = location_message->y;
			z2 = location_message->z;
			s2 = location_message->state;
			
			/* Calculate the distance between agent and message sending agent */
			/* Do not square root as more efficient to square interaction radius */
			distance = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
			/* If distance is equal or smaller than the agent interaction radius */
			if(distance < closestdist)
			{
				if(
				/* If MKK agent is free in cytoplasm and other agent is free cytoplasmic Trib or free cytoplasmic MAPK */
				(s1 == 0 && s2 == 900) ||
				(s1 == 0 && s2 == 200) ||

				/* if agent is free cytoplasmic MAPK and the other is reacted-Trib */
				(s1 == 2 && s2 == 403) ||
				(s1 == 2 && s2 == 400) ||				
				(s1 == 2 && s2 == 500) || 	
				/* if agent is free cytoplasmic MAPK and the other is reacted-Trib */
				(s1 == 200 && s2 == 1)
				)
				/* Find closest interaction partner */
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
		location_message = get_next_location_message(location_message);
		

	}
		
	/* If closest id is not the dummy value therefore a possible bond has been found */
	if(closestid != get_id())
	{
		//printf("%d(%d)\t> try bond message:\tto:%d(%d)\tdistance:%f\n", get_id(), get_state(), closestid, closeststate, closestdist);
		/* Send try to bond message */
		add_newbond_message(get_id(),         
						get_state(),
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
int Protein_checkbondtries()
{
	xmachine_memory_Protein * xmemory = current_xmachine->xmachine_Protein;
	
	double closestdist = 0.0;
	int closestid = xmemory->id; // dummy value
	int closeststate = 0;

	/* Check if I have been bound */
	/* Get first bond message */
	newbond_message = get_first_newbond_message();
	/* And loop through messages until there are none left */
	while(newbond_message)
	{
		if(newbond_message->idto == get_id() && newbond_message->bindunbind == 3)
		{
			if((closestid == xmemory->id) || (newbond_message->distance < closestdist))
			{
				closestdist = newbond_message->distance;
				closestid = newbond_message->idfrom;
				closeststate = newbond_message->statefrom;
			}
		}
		
		/* Get next trybond message */
		newbond_message = get_next_newbond_message(newbond_message);
	}
	// if not dummy value then bond created
	if(closestid != xmemory->id)
	{
		/* Send bond message to molecule to bind */
		add_finalbond_message(xmemory->id,
						xmemory->state,
						closestid,          /* id of agent message is refering to */
						0,                  /* 0 corresponds to a bind tag */
						0.0,    			/* distance not used */
						get_iradius(),
						get_posx(),
						get_posy(),
						get_posz());
		
		/* Then check state and change accordingly */
		/* If SB-TF bound to free pMK in nucleus change SB-TF to activated SB-TF */
		 if(get_state() == 400   && closeststate == 2)
			{
				set_state(21);
			}
		/* /* If SB-TF bound to free pMK in nucleus change MB-TF to activated MB-TF  */ 
		 if(get_state() == 403   && closeststate == 2)
			{
				set_state(22);
			}	
		/* If pMK free in nucleus and bound by exporting nuclear receptor then make pMK bound to exporting nuclear receptor [ExR], and change state to active MAPK */
		if(get_state() == 2  && closeststate == 500) 
			{
				set_state(200);						
			}
		/* If MAPK free in the cytoplasm and bound to MAPKK change state to pMK and transolcate into the nucleous */
		if(get_state() == 200   && closeststate == 0) 
			{
				set_state(2);
				set_mover(0.0); 
			}
		/* If free and active Trib in the cytoplasm and bound to MKK in cytoplasm change state to MT3 */
		if(get_state() == 900  && closeststate == 0)
			{
				set_state(1);
			}
				
	}	
		/* If free cytoplasmic reacted-Trib bind to free MAPK in the cytoplasm, change state to and re-set RADP4 to 50 */		
		if(get_state() == 1 && closeststate == 200)
			{
				set_state(11);
				set_RADP4(50);									
			}
	

	return 0;
}

int Protein_move()
{
	xmachine_memory_Protein * xmemory = current_xmachine->xmachine_Protein;
	
	double cellradius = 10000.0;    /* approximate cell radius in nano meters */
	double nuclradius = 4000.0;     /* approximate cell nucleus radius in nano meters */
	double dt = 10.0;                /* Time differential */
	/* Molecule movement constants */
	double speed_ave = 2;			/* Average speed of protein agent */
	double speed_range = 1;      /* Range of speed, + or - from average speed */
	double angle_range = pi/10.0;   /* Angle range in radians, + or - from current angle */
	/* Variables for calculations */
	double movex;
	double movey;
	double movez;

	/* Check if I have been bound */
	/* Get first bond message */
	finalbond_message = get_first_finalbond_message();
	/* And loop through messages until there are none left */
	

	while(finalbond_message)
	{
		/* bond message refers to me (same id) and tag refers to a bind */
		if(finalbond_message->idto == xmemory->id && finalbond_message->bindunbind == 0)
		{
			/* Then check state and change accordingly */
			/* If pMK free in nucleus and bound to SB-TF change pMK state to dormant-pMK bound to SB-TF  and re-set RADP to 100 */
			if(get_state() == 2   && finalbond_message->statefrom == 400) 
				{
					set_state(42); 
					set_RADP(100);
				}
			/* If pMK free in nucleus and bound to MB-TF change pMK state to dormant-gene expression activator (DGEA) pMK bound to MB-TF  and re-set RADP to 19 */
			if(get_state() == 2   && finalbond_message->statefrom == 403) 
				{
					set_state(24); 
					set_RADP(19);
				}
			/* free pMK binding to active Exrec change ExR state to dormant ExR*/
			if(get_state() == 2  && finalbond_message->statefrom == 500) set_state(501); 
			/* If MKK in cytoplasm and bound by Trib; change state to dormant MKK [dMKK] and re-set RADP to 7 */
			if(get_state() == 0   && finalbond_message->statefrom == 900) 
				{
					set_state(55); 
					set_RADP(7);
				}
			/* If free cytoplasmic MK then and bound to MT3; change MMT3 in cytoplasam and re-set RADP3 to 20 */
			if(get_state() == 200 && finalbond_message->statefrom == 1) 
				{
					set_state(222);
					set_RADP3(20);
				}
			/* If free cytoplasmic MKK then and bound to MK change state to dMKK free in the cytoplasm and re-set RADP to 7*/
			if(get_state() == 0 && finalbond_message->statefrom == 200) 
				{
					set_state(55);
					set_RADP(7);
				} 
								
		}
		/* Get next bond message */
		finalbond_message = get_next_finalbond_message(finalbond_message);
	}

	/* If freely moving molecule */
	if((xmemory->state==900         /* trib free in cyto */
			|| xmemory->state==0          /* MKK free in cytoplasm */
			|| xmemory->state== 1         /* MT3 free in cytoplasm */
			|| xmemory->state==200          /* MK free in cytoplasm */
			|| xmemory->state==11          /* MMT3 free in cytoplasm */
			|| xmemory->state==2          /* pMK free in nucleous*/
			|| xmemory->state==400			 /* dormant unbound AP-1 free in nucleous*/
			|| xmemory->state==401			/* dormant bound AP-1 free in nucleous*/
			|| xmemory->state==403			/* Active unbound AP-1 free in nucleous*/
			|| xmemory->state==404
			|| xmemory->state==55 			/* dormant pMKK free in cytoplasm */
			|| xmemory->state==24			 /* dormant pMK free in cytoplasm */
			|| xmemory->state==222 ))			/* dormant dMK free in cytoplasm*/
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
		
		/* Calculate redirection of molecules if leaving prescribed compartment */
		/* If in cytoplasm */
		if((xmemory->state == 0 			/* MKK free in cytoplasm */
			|| xmemory->state == 1 			/* MT3 free in cytoplasm */
			|| xmemory->state==900         /* trib free in cytoplasm */
			|| xmemory->state==200          /* MK free in cytoplasm */
			|| xmemory->state==11
			|| xmemory->state==5			/* dMKK free in cytoplasm */
			|| xmemory->state==222))          /*dMK free in cytoplasm */
		{
			/* If position is beyond the cell membrane */
			if(xmemory->posr > cellradius)
			{
				/* Mirror position back into cell */
				xmemory->posr = (2*cellradius)-xmemory->posr;
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
			/* If position is beyond the nuclear membrane */
			if(xmemory->posr < nuclradius)
			{
				/* Mirror position back into cytoplasm */
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
		/* If in nucleus */
		else
		{
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
	}
	
	return 0;
}

int Receptor_outputdata() /* AP-1 was treated like Imrec */
{
	/* Send location message */
	add_location_message(get_id(),		/* agent id */
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
				if(get_state() == 500) set_state(501);				
			}
		}
		/* If recdelay is zero, re-set redelay to 10*/
		if(get_recdelay() == 0) set_recdelay(10);

	
	return 0;
}
int Receptor_inputdata()
{
	xmachine_memory_Receptor * xmemory = current_xmachine->xmachine_Receptor;
	
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
	location_message = get_first_location_message();
	/* And loop through them until there are none left */
	while(location_message)
	{
		/* If the location message is not from the current agent */
		if((location_message->id != get_id()))
		{
			/* Copy location message data to #2 variables */
			id2 = location_message->id;
			x2 = location_message->x;
			y2 = location_message->y;
			z2 = location_message->z;
			s2 = location_message->state;
			
			/* Calculate the distance between agent and message sending agent */
			/* Do not square root as more efficient to square interaction radius */
			distance = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
			/* If agent is active Exporting Nuclear Receptor [ExR] */
			if(s1 == 500 && s2 == 2)
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
		/* Get next location message */
		location_message = get_next_location_message(location_message);
	}
		
	/* If closest id is not the dummy value therefore a possible bond has been found */
	if(closestid != xmemory->id)
	{
		/* Send try to bond message */
		add_newbond_message(xmemory->id,         
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
int Receptor_move()
{
	xmachine_memory_Receptor * xmemory = current_xmachine->xmachine_Receptor;
	
	double nuclradius = 4000.0;     /* Approximate cell nucleus radius in nano meters */		
	double dt = 10.0;                 /* Time differential */						
	/* Nuclear receptor constants */
	double nuc_speed_ave = 2.0;                /* Average speed of  ... */ 
	double nuc_speed_range = 1.0;               /* Range of speed, + or - from average speed */
	double nuc_ang_dev = pi/10.0;               /* Angle dev */  
	double nuc_ang_speed_ave = nuc_speed_ave;    /* Angle speed average */
	double nuc_ang_speed_dev = nuc_speed_range;  /* Angle speed dev */ 
	
	/* Check if I have been bound */
	/* Get first bond message */
	finalbond_message = get_first_finalbond_message();
	/* And loop through messages until there are none left */
	while(finalbond_message)
	{
		/* bond message refers to me (same id) and tag refers to a bind */
		if(finalbond_message->idto == xmemory->id && finalbond_message->bindunbind == 0)
		{
			/* If dorment exporting receptor */
			if(xmemory->state == 500)
			{
				/* Remember ID of molecule bound to me */
				xmemory->boundindex = finalbond_message->idfrom;
			}
			
			if(xmemory->state == 501)
			{
				/* Remember ID of molecule bound to me */
				xmemory->boundindex = finalbond_message->idfrom;
			}
		}
		/* Get next bond message */
		finalbond_message = get_next_finalbond_message(finalbond_message);
	}

	/* If nuclear receptor */
	
	if(xmemory->state > 401 && xmemory->state < 600)
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
