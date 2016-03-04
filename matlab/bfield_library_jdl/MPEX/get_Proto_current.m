function [helicon_current,current_A,current_B,config,skimmer] = get_Proto_current(shot)


switch shot
    case 7400
        helicon_current = 30; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
    case 7403
        helicon_current = 60; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
    case 7404
        helicon_current = 90; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
    case 7405
        helicon_current = 120; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
    case 7406
        helicon_current = 150; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
    case 7408
        helicon_current = 180; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
    case 7410
        helicon_current = 210; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
    case 7412
        helicon_current = 240; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
    case 7413
        helicon_current = 270; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
    case 7416
        helicon_current = 300; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
    case 7417
        helicon_current = 350; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
    case 7418
        helicon_current = 400; current_A = 6368; current_B = 0; config = 'flat'; skimmer = 0;
        
    case 7445
        helicon_current = 210; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 0;

    
    case 7477
        helicon_current = 30; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;
    case 7483
        helicon_current = 60; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;
    case 7487
        helicon_current = 90; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;
    case 7488
        helicon_current = 120; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;    
    case 7492
        helicon_current = 150; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;   
    case 7493
        helicon_current = 180; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;         
    case 7494
        helicon_current = 210; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;              
    case 7495
        helicon_current = 240; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;               
    case 7496
        helicon_current = 270; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;               
    case 7497
        helicon_current = 300; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;                        
    case 7498
        helicon_current = 350; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;              
    case 7500
        helicon_current = 400; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;                
    case 7501
        helicon_current = 500; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;                
    case 7503
        helicon_current = 600; current_A = 3300; current_B = 0; config = 'flat'; skimmer = 1;                                
    otherwise
        error(['Did not recognize shot: ',num2str(shot)])
end
if skimmer
    fprintf('This shot includes skimmer\n')
end
