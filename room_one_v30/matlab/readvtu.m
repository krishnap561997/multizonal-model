function [sln,lrp,tag,npart,time,step] = readvtu(infile,sln_num,lrp_num,tag_num)

%%  Read header parameters
    fid = fopen(infile,'r');
    str = textscan(fid,'%s',7,'Delimiter','\n');
    
    % Find time
    tim = str{1,1}{4,1};

    i = 74;
    ii = 0;
    flag = 0;
    while flag == 0
        i = i + 1;
        ii = ii + 1;

        if (tim(i) == '<')
            flag = 1;
        else
            aux(ii) = tim(i);
        end
    end
    time = str2double(aux);
    clear aux;

    % Find istep
    step = str{1,1}{5,1};
    i = 72;
    ii = 0;
    flag = 0;
    while flag == 0
        i = i + 1;
        ii = ii + 1;

        if (step(i) == '<')
            flag = 1;
        else
            aux(ii) = step(i);
        end
    end
    step = str2double(aux);
    clear aux;

    % Find total number of particles
    part = str{1,1}{7,1};
    i = 23;
    ii = 0;
    flag = 0;
    while flag == 0
        i = i + 1;
        ii = ii + 1;

        if (part(i) == '"')
            flag = 1;
        else
            aux(ii) = part(i);
        end
    end
    npart = str2double(aux);
    clear aux;
    fclose(fid);


%%  Read solution data

    fprintf('----------------------------- \n');
    fprintf('Reading file: %s \n',infile);
    fprintf('Time = %f \n',time);
    fprintf('Step = %i \n',step);
%     fprintf(' \n');

    fid = fopen(infile,'r');

    % Find start of solution arrays
    flag=0;
    i = 0;
    while flag==0
        i = i + 1;
        A = fread(fid,1,'char*1');
        header(i) = char(A');
        if (header(i) == '_') && (header(i-2) == '>')
            flag = 1;
        end
    end

%    
    pos_len = fread(fid,1,'integer*4');
    pos = fread(fid,pos_len/4,'float32');

%
    sln_len = zeros(sln_num,1);
    sln = zeros(npart,sln_num);
    for i=1:sln_num
        sln_len(i) = fread(fid,1,'integer*4');
        sln(:,i) = fread(fid,sln_len(i)/4,'float32');
    end

%
    lrp_len = zeros(lrp_num,1);
    lrp = zeros(npart,lrp_num);
    for i=1:lrp_num
        lrp_len(i) = fread(fid,1,'integer*4');
        lrp(:,i) = fread(fid,lrp_len(i)/4,'float32');

    end

%
    tag_len = zeros(tag_num,1);
    tag = zeros(npart,tag_num);
    for i=1:tag_num
        tag_len(i) = fread(fid,1,'integer*4');
        tag(:,i) = fread(fid,tag_len(i)/4,'float32');
    end
    
    fclose(fid);

end
