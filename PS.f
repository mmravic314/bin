c
c      Fit a great circle arc to a set of points
c
c      Brian J. Smith
c
c      Department of Chemistry, La Trobe Institute for Molecular Science,
c      La Trobe University, Melbourne, Victoria 3086, Australia.
c      &
c      The Walter & Eliza Hall Institute of Medical Research,
c      Parkville, Victoria 3052, Australia.
c
c      Use of any part of this program should cite
c      B.J. Smith, PS – a program for the analysis of helix geometry, 
c      J. Mol. Graph. Model. (2011), doi:10.1016/j.jmgm.2011.11.004
c
c      No warranties, implied or otherwise, accompany this program.
c
c      Usage: PS -i filename -o output -H helix -c chain -f first 
c                -l last -a atom_type -m model -r ref_type 
c                -k conformer -R Radius -s seed -S scan -h 
c
c      filename  = name of input pdb file (required)
c      output    = name of output file
c      helix     = index of helix in pdb file
c      chain     = chain identifier
c      first     = index of first residue
c      last      = index of last residue
c      model     = model number in multi-structure NMR PDB file
c      atom_type = name of atom(s) to select (up to five)
c      ref_type  = name of atom(s) for comparison (up to five)
c      conformer = choice of conformer 
c      Radius    = fixed Radius (Angstrom)
c      seed      = initial radius (Angstrom)
c      scan      = fine/medium/coarse
c
c      Maximum number of points = 500
c      Maximum number of atom types = 5
c
       real*8        x(500),y(500),z(500)
       real*8        u(500),v(500),w(500)
       real*8        dist(500)
       real*8        objective,param(7),grad(7)
       real*8        dref(500)
       real*8        rise(500),avrise,varrise,sdrise
       real*8        theta(500),avtheta,vartheta,sdtheta
       real*8        xr,yr,zr
       real*8        d,dav,maxdev,mindev,s2,sd,Radius,seed
       real*8        length,pitch,Len,norm,step
       real*8        pi
       real*8        best_obj,best_seed
       real*8        dir(6)
       integer*4     resnum(500),resndx(500),rnum(500)
       integer*4     arg,nargs
       integer*4     first,last,model,ModNum,pdbhelix
       integer*4     i,j,k,l
       integer*4     natoms,resid,natom_types,nref_types,nrefs,nconfs
       integer*4     ier,out
       integer*4     imax,imin,nit,MAXATM
       integer*4     fix_radius
       integer*4     scan_interval
       character*132 argument
       character*132 filename,output
       character*80  String
       character*80  atom_type,ref_type
       character*80  first_res,last_res
       character*12  rname(500)
       character*6   scan_type
       character*4   atname(500)
       character*4   atoms(5),refs(5)
       character*3   resname(500)
       character*1   conf(500),chid(500),seqins(500)
       character*1   chain
       character*1   first_ins,last_ins
       character*1   char,konform
       logical*1     default_chain
c
       common/coord/x,y,z,u,v,w,natoms
       common/cg/objective,param,grad,fix_radius
       common/io/out
c
c      Initialise some variables
c
       default_chain=.false.
       chain(1:1)=' '
       konform(1:1)=' '
       first_ins(1:1)=' '
       last_ins(1:1)=' '
       scan_type(1:6)='      '
       do i=1,80
        atom_type(i:i)=' '
        ref_type(i:i)=' '
       enddo
       ier=999  !  Error message from conjugate gradient minimizer
       out=6  !  Output unit
       pdbhelix=0
       nconfs=1
       first=99999
       last=-99999
       model=-99
       ModNum=-99
       nref_types=0
       fix_radius=1
       MAXATM=500
       Radius=0.0D0
       seed=0.0D0  !  Estimate seed
       best_seed=seed
       best_obj=99999.99999D10
       pi=4.0D0*atan(1.0D0)
c
       arg=0
       nargs=26
c
c      Parse input string
c
       do while (arg.le.nargs)
       call getarg(arg,argument)
        if (argument.eq.'-i') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(A)') filename
         open(unit=4,file=filename,status='old')
        elseif (argument.eq.'-H') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(I3)') pdbhelix
        elseif (argument.eq.'-o') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(A)') output
         open(unit=7,file=output,status='new')
         out=7
        elseif (argument.eq.'-c') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(A)') chain
        elseif (argument.eq.'-f') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(A80)') first_res
         call rhj(first_res)
         call res2res(first_res,first,first_ins)
        elseif (argument.eq.'-l') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(A80)') last_res
         call rhj(last_res)
         call res2res(last_res,last,last_ins)
        elseif (argument.eq.'-a') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(A80)') atom_type
        elseif (argument.eq.'-r') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(A80)') ref_type
         nref_types=1
        elseif (argument.eq.'-m') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(I3)') model
        elseif (argument.eq.'-k') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(A)') konform
        elseif (argument.eq.'-R') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(1F10.3)') Radius
         fix_radius=0
        elseif (argument.eq.'-s') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(1F10.3)') seed
        elseif (argument.eq.'-S') then
         arg=arg+1
         call getarg(arg,argument)
         read(argument,'(A)') scan_type
        elseif ((argument.eq.'-h').or.
     +          ((argument.ne.'-i').and.
     +           (argument.ne.'-o').and.
     +           (argument.ne.'-H').and.
     +           (argument.ne.'-c').and.
     +           (argument.ne.'-f').and.
     +           (argument.ne.'-l').and.
     +           (argument.ne.'-a').and.
     +           (argument.ne.'-m').and.
     +           (argument.ne.'-k').and.
     +           (argument.ne.'-r').and.
     +           (argument.ne.'-R').and.
     +           (argument.ne.'-s').and.
     +           (arg.ne.0).and.
     +           (argument(1:1).ne.' '))) then
       write(6,*)arg,argument
         write(6,'(A)')"Usage: PS -i filename -o output -H helix"
         write(6,'(10X,A)')"-c chain -f first -l last -a atom_type"
         write(6,'(10X,A)')"-m model -k conformer -r ref_type"
         write(6,'(10X,A)')"-R radius -s seed -S scan -h"
         write(6,'(A)')" "
         write(6,'(A)')" filename  = name of pdb file (A)"
         write(6,'(A)')" output    = name of output file (A)"
         write(6,'(A)')" helix     = index of helix in pdb file (I)"
         write(6,'(A)')" chain     = chain identifier (A)"
         write(6,'(A)')" first     = index of first residue (A)"
         write(6,'(A)')" last      = index of last residue (A)"
         write(6,'(A)')" model     = model number (in multi-structure NM
     +R file) (I)"
         write(6,'(A)')" conformer = choice of conformer - default is al
     +l (A)"
         write(6,'(A)')" atom_type = name of atom(s) to select (comma se
     +parated) (A)"
         write(6,'(A)')" ref_type  = name of atom(s) for comparison (com
     +ma separated) (A)"
         write(6,'(A)')" Radius    = fixed Radius (Angstrom) (F)"
         write(6,'(A)')" seed      = initial radius for CG optimiser (An
     +gstrom) (F)"
         write(6,'(A)')" scan      = scan type (fine/medium/coarse) (A)"
         write(6,'(A)')" "
         write(6,'(A)')" (A) - string, (I) - integer, (F) - float"
         call exit(0)
         stop
        endif
        arg=arg+1
       enddo
c
c      Retrieve information from HELIX record in pdb file
c
       if(pdbhelix.ne.0) then
        do while (.true.)
         read(4,'(A)',end=1)String
         if(String(1:5).eq.'HELIX') then
          read(String(8:10),'(I3)')i
          if(i.eq.pdbhelix) then
           read(String(20:20),'(A)')chain
           read(String(22:25),'(I4)')first
           read(String(26:26),'(A)')first_ins
           read(String(34:37),'(I4)')last
           read(String(38:38),'(A)')last_ins
           goto 2
          endif
         endif
        enddo
c
c       HELIX card not found
c
 1      continue
        write(6,'(A)')"HELIX card not found."
        write(6,'(A)')"Program must stop"
        call exit(1)
        stop
 2      continue
        rewind(4)
       endif
c
c      Check the consistency of the input
c
       if(first.gt.last) then
        write(out,'(A)')"Index of first (-f) residue is greater than las
     +t (-l)"
        write(out,'(A8,I5,A1,A8,I5,A1)')"first = ",first,first_ins,
     +                                  " last = ",last,last_ins
        write(out,'(A)')"Program must stop"
        call exit(2)
        stop
       endif
c
c      Check for a valid scan interval
c
       if((scan_type(1:6).ne.'      ').and.
     +    (scan_type(1:4).ne.'fine').and.
     +    (scan_type(1:6).ne.'medium').and.
     +    (scan_type(1:6).ne.'coarse')) then
        write(out,'(A)')"Unable to read scan type"
        write(out,'(A)')"Program must stop"
        call exit(4)
        stop
       endif
c
c      Write an acknowledgement
c
       write(out,'(A)')" PS (pandus semita) - a helix analysis program. 
     +Vers 1.0"
       write(out,'(A)')" "
       call get_command(argument)
       write(out,'(A)')"Invoked using: ",argument
c
c      Apply default atom type if necessary
c
       if(atom_type(1:4).eq.'    ') then
        write(out,'(A)')" "
        write(out,'(A)')" No atom type supplied - assuming CA."
        atom_type(1:4)='CA  '
       endif
c
c      Unpack the string of (comma separated) atom types (up to 5 types)
c
       natom_types=1
       atoms(1)(1:4)='    '
       i=0
       do j=1,80
        read(atom_type(j:j),'(A1)')char
        if(char(1:1).eq.',') then
         natom_types=natom_types+1
         atoms(natom_types)(1:4)='    '
         i=0
        elseif(char(1:1).eq.' ') then
         goto 3
        else
         i=i+1
         atoms(natom_types)(i:i)=char(1:1)
        endif
       enddo
 3    continue
c
c      Unpack the string of (comma separated) reference atom types (up to 5 types)
c
      if(nref_types.gt.0) then
       nref_types=1
       refs(1)(1:4)='    '
       i=0
       do j=1,80
        read(ref_type(j:j),'(A1)')char
        if(char(1:1).eq.',') then
         nref_types=nref_types+1
         refs(nref_types)(1:4)='    '
         i=0
        elseif(char(1:1).eq.' ') then
         goto 4
        else
         i=i+1
         refs(nref_types)(i:i)=char(1:1)
        endif
       enddo
      endif
 4    continue
c
       if(pdbhelix.ne.0) then
        write(out,'(A)')" "
        write(out,'(A)')" HELIX record read."
       endif
c
c      Extract data from the PDB file
c
       natoms=0
       do while (.true.)
        read(4,'(A)',end=5)String
c
        if(String(1:5).eq.'MODEL') then
          if(model.eq.-99) then
           write(out,'(A)')" "
           write(out,'(A)')"Multi-structure NMR file - no model selected
     +. Assuming model 1."
          model=1
          endif
          read(String(6:14),'(I10)')ModNum
        endif
c
        if(((String(1:4).eq.'ATOM').or.(String(1:6).eq.'HETATM')).and.
     +     (model.eq.ModNum)) then
c
c      The first chain is the default if nothing else is specified
c
         if((chain(1:1).eq.' ').and.(.not.default_chain)) then
          chain(1:1)=String(22:22)
          write(out,'(A)')" "
          write(out,'(A29,A1,A2)')" No chain specified - using >",
     +    chain(1:1),"<."
          default_chain=.true.
         endif
         if(String(22:22).eq.chain(1:1)) then
          read(String(23:26),'(I4)') resid
           if( ((resid.gt.first).or.
     +          ((resid.eq.first).and.
     +           (first_ins(1:1).eq.' ').and.
     +           (String(27:27).eq.' ')).or.
     +          ((resid.eq.first).and.
     +           (first_ins(1:1).ne.' ').and.
     +           (ichar(String(27:27)).ge.ichar(first_ins)))) 
     +         .and.
     +         ((resid.lt.last).or.
     +          ((resid.eq.last).and.
     +           (last_ins(1:1).eq.' ').and.
     +           (String(27:27).eq.' ')).or.
     +          ((resid.eq.last).and.
     +           (last_ins(1:1).ne.' ').and.
     +           (ichar(String(27:27)).le.ichar(last_ins))))
     +         .and.
     +         ((konform(1:1).eq.' ').or.
     +          ((konform(1:1).ne.' ').and.
     +          ((String(17:17).eq.konform(1:1)).or.
     +           (String(17:17).eq.' ')))) ) then
c
           call lhj(String(13:16))
           do j=1,natom_types
            if(String(13:16).eq.atoms(j)(1:4)) then
             natoms=natoms+1
c
c     Check array boundary
c
             if(natoms.ge.MAXATM) then
              write(out,'(A11,I4,A15)')"Maximum of ",MAXATM,
     +         " atoms exceeded"
              write(out,'(A)')"Program must stop"
              call exit(3)
              stop
             endif
c
c     Helix parameters (rise, phase yield) are calculated for
c     the first atom type, and first conformer ONLY. The first
c     conformer is designated 'A'.
c
             if(String(17:17).ne.' ') then
              nconfs=2 ! i.e. more than one conformation
             endif
             resndx(natoms)=0
             if(j.eq.1) then
              if((String(17:17).eq.' ').or.
     +           (String(17:17).eq.konform(1:1)).or.
     +          ((String(17:17).eq.'A').and.
     +           (konform(1:1).eq.' '))) then
               resndx(natoms)=1
              endif
             endif
c
             read(String(13:16),'(A4)')atname(natoms)
             read(String(17:17),'(A1)')conf(natoms)
             read(String(18:20),'(A3)')resname(natoms)
             read(String(22:22),'(A1)')chid(natoms)
             read(String(23:26),'(I4)')resnum(natoms)
             read(String(27:27),'(A1)')seqins(natoms)
             read(String(31:54),'(3F8.3)')
     +        x(natoms),y(natoms),z(natoms)
            endif
           enddo
          endif
         endif
        endif
       enddo
 5     continue
       if(((fix_radius.eq.1).and.(natoms.lt.7)).or.
     +    ((fix_radius.eq.0).and.(natoms.lt.6)))then
        write(out,'(A)')"System overdetermined - need more than 6 points
     + (or 5 points with fixed Radius) in CG minimization"
        write(out,'(A)')"Program must stop"
        call exit(5)
        stop
       endif
       if(natoms.eq.0) then
        write(out,'(A)')"No atoms read"
        if((model.ne.-99).and.(ModNum.eq.-99)) then
         write(out,'(A)')"Model (-m) selected but no models found in PDB
     + file"
        endif
        write(out,'(A)')"Program must stop"
        call exit(6)
        stop
       endif
c
c      Scan the seed to find the optimum
c
       if(scan_type(1:1).ne.' ') then
        if(scan_type(1:4).eq.'fine') then
         scan_interval=20
        elseif(scan_type(1:6).eq.'medium') then
         scan_interval=98
        elseif(scan_type(1:6).eq.'coarse') then
         scan_interval=196
        else
         write(out,'(A)')"Unable to read scan interval"
         write(out,'(A)')"Program must stop"
         call exit(7)
         stop
        endif 
        write(out,'(A)')" "
        write(out,'(A17,A6,A2,$)')
     +   " Scanning seed - ",scan_type(1:6),": "
        write(out,'(A)')" "
        do i=20,1000,scan_interval
         ier=999
c
c      Initialise the parameters
c
         param(7)=i
         call initialize
c
         call acgm(0,1.0D-4,ier,nit,step)
         call acgm(1,1.0D-12,ier,nit,step)
         if(objective.le.best_obj) then
          best_obj=objective
          best_seed=i
          write(out,'(A,$)')"*"
         else
          write(out,'(A,$)')"."
         endif
        enddo
        write(out,'(A)')" "
       endif
c
c      Repeat with best seed if one found in scan, 
c      else use the user specified seed or estimate
c
       if(fix_radius.eq.1) then
        param(7)=best_seed
       else
        param(7)=Radius
       endif
       call initialize
c
c      Conjugate gradient minimization (of sum dist**2)
c
       call acgm(0,1.0D-4,ier,nit,step)
c
c      Conjugate gradients minimization (of variance)
c
       call acgm(1,1.0D-12,ier,nit,step)
c
c      Check for sensible sphere radius
c
        if(param(7).le.0.0D0) then
         write(out,'(A21,1F8.3,A1)')
     +    "Radius not positive (",param(7),")"
         write(out,'(A)')"Program must stop"
         call exit(8)
         stop
        endif
c
c      Report the optimized axis parameters if mimimization successful
c      Print a warning if termination not successful
c
       if(ier.ne.0) then
        write(out,'(A)')"Warning: CG optimisation did not terminate succ
     +esfully"
        write(out,'(A)')"The results provided below should be scrutinise
     +d"
        write(out,'(A)')"Pay particular attention to the residual gradie
     +nts"
       endif
        write(out,'(A)')" "
        write(out,'(A,1E12.6,/)')" Minimized Objective : ",objective
        write(out,'(A)')
     +  " Final parameters (residual gradients in parenthesis)"
        write(out,'(25X,A1)')"R"
        write(out,'(A16,1F12.3,A2,1E8.1,A1)')
     +  " Sphere radius: ",param(7)," (",grad(7),")"
        write(out,'(25X,A1,20X,A1,20X,A1)')"a","b","c"
        write(out,'(A20,3(F8.3,A2,1E8.1,A3))')" Sphere centre:     ",
     +  param(1)," (",grad(1),")  ",
     +  param(2)," (",grad(2),")  ",
     +  param(3)," (",grad(3),")  "
        write(out,'(25X,A1,20X,A1,20X,A1)')"l","m","n"
c
c       Normalize the plane vectors
c
        norm=dsqrt(param(4)**2 + param(5)**2 + param(6)**2)
        param(4)=param(4)/norm
        param(5)=param(5)/norm
        param(6)=param(6)/norm
c
        write(out,'(A20,3(F8.3,A2,1E8.1,A3))')" Plane orient:      ",
     +  param(4)," (",grad(4),")  ",
     +  param(5)," (",grad(5),")  ",
     +  param(6)," (",grad(6),")  "
        write(out,'(A)')" "
c
c      Determine location of axis perpendiculars
c
        do i=1,natoms
         call UVW(param(1),param(2),param(3),x(i),y(i),z(i),
     +     param(4),param(5),param(6),param(7),u(i),v(i),w(i))
        enddo
c
c      Determine helix direction at terminii
c
        call helixdir(dir)
c
c      Report distances of points to helix axis
c
        dav=0.0D0
        maxdev=-99999.99999D0
        mindev=99999.99999D0
        imax=-9999
        imin=9999
        write(out,'(34X,A27)')"X       Y       Z      Dist"
        do i=1,natoms
         call arcdist(param(1),param(2),param(3),dist(i),
     +   x(i),y(i),z(i),param(4),param(5),param(6),param(7))
         write(out,'(8X,I3,1X,A4,A1,A3,1X,A1,I4,A1,3X,4F8.3)')
     +   i,atname(i),conf(i),resname(i),chid(i),resnum(i),
     +   seqins(i),x(i),y(i),z(i),dist(i)
         dav=dav+dist(i)
         if(dist(i).ge.maxdev) then
          maxdev=dist(i)
          imax=i
         endif
         if(dist(i).le.mindev) then
          mindev=dist(i)
          imin=i
         endif
        enddo
        dav=dav/float(natoms)
        s2=0.0D0
        do i=1,natoms
         s2=s2+(dist(i)-dav)**2
        enddo
        s2=dsqrt(s2/(float(natoms-1)))
        write(out,'(A)')" "
        write(out,'(A20,1F8.3,A5,1F5.3,A9)')" Average distance = ",
     +  dav," +/- ",s2,' Angstrom'
        write(out,'(A20,1F8.3,A19,A3,I4)')" Maximum distance = ",
     +  maxdev,"           Residue ",resname(imax),imax+first-1
        write(out,'(A20,1F8.3,A19,A3,I4)')" Minimum distance = ",
     +  mindev,"           Residue ",resname(imin),imin+first-1
c
c     Report locations of axis perpendiculars
c
        write(out,'(A)')" "
        write(out,'(A)')"Axis perpendiculars:"
        write(out,'(64X,A9)')"Direction"
        write(out,'(34X,A45)')"X       Y       Z       alpha   beta    g
     +amma"
        do i=1,natoms
         if(i.eq.1) then
          write(out,'(A6,2X,I3,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,1X,
     +                3F8.3)')
     +     "HETATM",i,"X   ",conf(i),resname(i),chid(i),resnum(i),
     +      seqins(i),u(i),v(i),w(i),dir(1),dir(2),dir(3)
         elseif(i.eq.natoms) then
          write(out,'(A6,2X,I3,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,1X,
     +                3F8.3)')
     +     "HETATM",i,"X   ",conf(i),resname(i),chid(i),resnum(i),
     +      seqins(i),u(i),v(i),w(i),dir(4),dir(5),dir(6)
         else
          write(out,'(A6,2X,I3,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3)')
     +     "HETATM",i,"X   ",conf(i),resname(i),chid(i),resnum(i),
     +      seqins(i),u(i),v(i),w(i)
         endif
        enddo
c
c     Collect only interesting atoms
c
      j=0
      k=0
      do i=1,natoms
       if(seqins(i)(1:1).eq.' ') then
        seqins(i)='@'
       endif
       if(resndx(i).ne.0) then
        j=j+1
        if(j.eq.1) then
         rnum(j)=1
        else
         rnum(j)=rnum(j-1)+
     +    (ichar(seqins(i))-ichar('@'))+
     +    (resnum(i)-resnum(k))
        endif
        x(j)=x(i)
        y(j)=y(i)
        z(j)=z(i)
        u(j)=u(i)
        v(j)=v(i)
        w(j)=w(i)
        dist(j)=dist(i)
c
        rname(j)(1:1)=conf(i)(1:1)
        rname(j)(2:4)=resname(i)(1:3)
        if(seqins(i)(1:1).eq.'@') then
         seqins(i)(1:1)=' '
        endif
        rname(j)(5:5)=' '
        rname(j)(6:6)=chid(i)(1:1)
        rname(j)(7:7)=' '
        write(rname(j)(8:11),'(I4)')resnum(i)
        rname(j)(12:12)=seqins(i)(1:1)
c
        k=i
       endif
      enddo
      natoms=j
c
c     Determine rise per residue pair
c
       avrise=0.0D0
       do i=1,natoms-1
         rise(i)=dsqrt((u(i+1)-u(i))**2+
     +                 (v(i+1)-v(i))**2+
     +                 (w(i+1)-w(i))**2)
         rise(i)=rise(i)/(rnum(i+1)-rnum(i))
         avrise=avrise+rise(i)
       enddo
       avrise=avrise/float(natoms-1)
       varrise=0.0D0
       do i=1,natoms-1
        varrise=varrise+((avrise-rise(i))**2)
       enddo
       varrise=varrise/float(natoms-2)
       sdrise=dsqrt(varrise)
c
c     Determine helical phase yield per residue pair
c     
       avtheta=0.0D0
       do i=1,natoms-1
         costheta=((x(i)-u(i))*(x(i+1)-u(i+1))+
     +             (y(i)-v(i))*(y(i+1)-v(i+1))+
     +             (z(i)-w(i))*(z(i+1)-w(i+1)))/
     +      (dsqrt((x(i)-u(i))**2+
     +             (y(i)-v(i))**2+
     +             (z(i)-w(i))**2)*
     +       dsqrt((x(i+1)-u(i+1))**2+
     +             (y(i+1)-v(i+1))**2+
     +             (z(i+1)-w(i+1))**2))
         theta(i)=acos(costheta)/(rnum(i+1)-rnum(i))
         avtheta=avtheta+theta(i)
       enddo
       avtheta=avtheta/float(natoms-1)
       vartheta=0.0D0
       do i=1,natoms-1
        vartheta=vartheta+(avtheta-theta(i))**2
       enddo
       vartheta=vartheta/float(natoms-2)
       sdtheta=dsqrt(vartheta)
c
       if((natom_types.gt.1).or.(nconfs.gt.1)) then
        write(out,'(A)')" "
        write(out,'(A)')"Helix parameters (rise, phase yield)"
        write(out,'(A)')"are calculated for the first atom type and conf
     +ormer ONLY."
       endif
c
       write(out,'(A)')" "
       write(out,'(30X,A)')"  Residue     Rise     Phase"
       write(out,'(30X,A)')" separation            yield"
       write(out,'(30X,A)')"           (Angstrom)(degrees)"
       do i=1,natoms-1
        write(out,'(I3,2X,A12,A4,A12,2X,I1,4X,1F8.3,3X,1F8.3)')
     +  i,rname(i)," -> ",rname(i+1),rnum(i+1)-rnum(i),rise(i),
     +  180.0D0*theta(i)/pi
       enddo
       write(out,'(A)')" "
       write(out,'(A,1F8.3,A,1F6.3,A9)')
     +  ' Rise/residue: ',avrise,
     +  ' +/- ',sdrise,' Angstrom'
       write(out,'(A15,1F8.3,A5,1F6.3,A8)')
     +  ' Phase:        ',180.0D0*avtheta/pi,
     +  ' +/- ',(180.0D0*sdtheta/pi),
     +  ' degrees'
       write(out,'(A15,1F8.3,A5,1F6.3)')
     +  ' Residue/turn: ',2.00D0*pi/avtheta,
     +  ' +/- ',2.0D0*pi*sdtheta/(avtheta**2)
c
c      Determine helix length
c
       length=2.0D0*param(7)*
     +  asin(dsqrt((u(natoms)-u(1))**2+
     +             (v(natoms)-v(1))**2+
     +             (w(natoms)-w(1))**2)/(2.0D0*param(7)))
c
c      Determine helix pitch
c
       pitch=2.0D0*pi*length/(float(natoms-1)*avtheta)
       write(out,'(A15,1F8.3,A5,1F6.3,A9)')
     +  ' Pitch:        ',pitch,
     +  ' +/- ',
     +  2.0D0*pi*length*dsqrt(vartheta)/(float(natoms-1)*avtheta**2),
     +  ' Angstrom'
c
       write(out,'(A15,1F8.3,A20)')" Helix length: ",length,
     +  "            Angstrom"
c
c      Analyse the reference data
c
       if(nref_types.gt.0) then 
        do j=1,nref_types 
         write(out,'(A)')" "
         write(out,'(A30,A4)')"Reference data for atom type: ",refs(j)
          write(out,'(26X,A43)')
     +    "X       Y       Z     Dist    Phase    Rise"
          write(out,'(26X,A43)')
     +    "                              yield        "
c
c      Initialize
c 
         dav=0.0D0
         maxdev=-99999.99999
         mindev=99999.99999
         imax=-9999
         imin=-9999
         nrefs=0
         i=0
         nconfs=1
c
         rewind(4)
         do while (.true.) 
          read(4,'(A)',end=6)String
c
          if(String(1:5).eq.'MODEL') then 
            read(String(6:14),'(I10)')ModNum
          endif 
c
          if(((String(1:4).eq.'ATOM').or.(String(1:6).eq.'HETATM')).and.
     +       (model.eq.ModNum)) then
           if(String(22:22).eq.chain(1:1)) then 
            read(String(23:26),'(I4)') resid
            if( ((resid.gt.first).or. 
     +           ((resid.eq.first).and.
     +            (first_ins(1:1).eq.' ').and.
     +            (String(27:27).eq.' ')).or.
     +           ((resid.eq.first).and.
     +            (first_ins(1:1).ne.' ').and.
     +            (ichar(String(27:27)).ge.ichar(first_ins))))
     +          .and.
     +          ((resid.lt.last).or.
     +           ((resid.eq.last).and.
     +            (last_ins(1:1).eq.' ').and.
     +            (String(27:27).eq.' ')).or.
     +           ((resid.eq.last).and.
     +            (last_ins(1:1).ne.' ').and.
     +            (ichar(String(27:27)).le.ichar(last_ins))))
     +          .and.
     +          ( (konform(1:1).eq.' ').or.
     +           ((konform(1:1).ne.' ').and.
     +           ((String(17:17).eq.konform(1:1)).or.
     +            (String(17:17).eq.' ')))) ) then
c
              call lhj(String(13:16))
              if(String(13:16).eq.refs(j)(1:4)) then 
               i=i+1
c
               if(String(17:17).ne.' ') then 
                nconfs=2 ! more than one
               endif 
               resndx(i)=0
                if((String(17:17).eq.' ').or. 
     +             (String(17:17).eq.konform(1:1)).or.
     +            ((String(17:17).eq.'A').and.
     +             (konform(1:1).eq.' ')))then
                 resndx(i)=1
                endif 
c
               read(String(13:16),'(A4)')atname(i)
               read(String(17:17),'(A1)')conf(i)
               read(String(18:20),'(A3)')resname(i)
               read(String(22:22),'(A1)')chid(i)
               read(String(23:26),'(I4)')resnum(i)
               read(String(27:27),'(A1)')seqins(i)
               read(String(31:54),'(3F8.3)')x(i),y(i),z(i)
               call arcdist(param(1),param(2),param(3),d,
     +          x(i),y(i),z(i),param(4),param(5),param(6),
     +          param(7))
               dref(i)=d
c - perpendiculars
               call UVW(param(1),param(2),param(3),
     +          x(i),y(i),z(i),
     +          param(4),param(5),param(6),param(7),
     +          u(i),v(i),w(i))
c
               dav=dav+d
               if(d.ge.maxdev) then 
                maxdev=d
                imax=resid
               endif 
               if(d.le.mindev) then 
                mindev=d
                imin=resid
               endif 
c
               nrefs=i
c
             endif 
            endif 
           endif 
          endif 
         enddo 
c
 6      continue
c
c - theta and rise
c
c     Collect only interesting atoms
c
         l=0
         k=0
         do i=1,nrefs 
          if(seqins(i)(1:1).eq.' ') then 
           seqins(i)='@'
          endif 
          if(resndx(i).ne.0) then 
           l=l+1
           if(l.eq.1) then 
            rnum(l)=1
           else 
            rnum(l)=rnum(l-1)+
     +       (ichar(seqins(i))-ichar('@'))+
     +       (resnum(i)-resnum(k))
           endif 
           atname(l)=atname(i)
           conf(l)=conf(i)
           resname(l)=resname(i)
           chid(l)=chid(i)
           resnum(l)=resnum(i)
           seqins(l)=seqins(i)
           x(l)=x(i)
           y(l)=y(i)
           z(l)=z(i)
           u(l)=u(i)
           v(l)=v(i)
           w(l)=w(i)
           dist(l)=dist(i)
c
          if(seqins(i)(1:1).eq.'@') then 
           seqins(i)(1:1)=' '
          endif 
c
           k=i
          endif 
         enddo 
         nrefs=l
c
         avtheta=0.0D0
         avrise=0.0D0
         do i=1,nrefs-1 
          costheta=((x(i)-u(i))*(x(i+1)-u(i+1))+
     +              (y(i)-v(i))*(y(i+1)-v(i+1))+
     +              (z(i)-w(i))*(z(i+1)-w(i+1)))/
     +       (dsqrt((x(i)-u(i))**2+
     +              (y(i)-v(i))**2+
     +              (z(i)-w(i))**2)*
     +        dsqrt((x(i+1)-u(i+1))**2+
     +              (y(i+1)-v(i+1))**2+
     +              (z(i+1)-w(i+1))**2))
          theta(i)=acos(costheta)/(rnum(i+1)-rnum(i))
          avtheta=avtheta+theta(i)
          rise(i)=dsqrt((u(i+1)-u(i))**2+
     +                  (v(i+1)-v(i))**2+
     +                  (w(i+1)-w(i))**2)
          rise(i)=rise(i)/(rnum(i+1)-rnum(i))
          avrise=avrise+rise(i)
         enddo 
         avtheta=avtheta/float(nrefs-1)
         avrise=avrise/float(nrefs-1)
         do i=1,nrefs-1 
          write(out,
     +    '(I3,1X,A4,1X,A1,A3,1X,A1,I4,A1,1X,4F8.3,1X,3F8.3)')
     +     i,atname(i),conf(i),resname(i),chid(i),resnum(i),seqins(i),
     +     x(i),y(i),z(i),dref(i),180.0D0*theta(i)/pi,rise(i)
         enddo 
         write(out,
     +    '(I3,1X,A4,1X,A1,A3,1X,A1,I4,A1,1X,4F8.3)')
     +    nrefs,atname(nrefs),conf(nrefs),resname(nrefs),chid(nrefs),
     +    resnum(nrefs),seqins(nrefs),x(nrefs),y(nrefs),z(nrefs),
     +    dref(nrefs)
c
         dav=dav/float(nrefs)
         s2=0.0D0
         vartheta=0.0D0
         varrise=0.0D0
         do i=1,nrefs
          s2=s2+(dref(i)-dav)**2
         enddo
         s2=s2/(float(nrefs-1))
         sd=dsqrt(s2)
         do i=1,nrefs-1 
          vartheta=vartheta+(avtheta-theta(i))**2
          varrise=varrise+(avrise-rise(i))**2
         enddo 
         vartheta=vartheta/float(nrefs-2)
         varrise=varrise/float(nrefs-2)
         sdtheta=dsqrt(vartheta)
         sdrise=dsqrt(varrise)
c
c      Determine helix length
c
         length=2.0D0*param(7)*
     +    asin(dsqrt((u(nrefs)-u(1))**2+
     +               (v(nrefs)-v(1))**2+
     +               (w(nrefs)-w(1))**2)/(2.0D0*param(7)))
c
c      Determine helix pitch
c
         pitch=2.0D0*pi*length/(float(nrefs-1)*avtheta)
c
         write(out,'(A)')" "
         write(out,'(A22,1F8.3,A5,1F6.3,A9)')" Average distance   = ",
     +   dav," +/- ",sd,' Angstrom'
         write(out,'(A22,1F8.3,A20,A3,I4)')" Maximum distance   = ",
     +   maxdev,"            Residue ",resname(imax-first+1),imax
         write(out,'(A22,1F8.3,A20,A3,I4)')" Minimum distance   = ",
     +   mindev,"            Residue ",resname(imin-first+1),imin
         write(out,'(A22,1F8.3,A5,1F6.3,A8)')" Average phase      = ",
     +   180.0D0*avtheta/pi," +/- ",(180.0D0*sdtheta/pi),' degrees'
         write(out,'(A22,1F8.3,A5,1F6.3,A9)')" Average rise       = ",
     +   avrise," +/- ",sdrise,' Angstrom'
         write(out,'(A22,1F8.3,A5,1F6.3,A9)')" Average pitch      = ",
     +   pitch," +/- ",
     +   2.0D0*pi*length*dsqrt(vartheta)/(float(nrefs-1)*avtheta**2),
     +   ' Angstrom'
         write(out,'(A22,1F8.3,A5,1F6.3)')" Av. residues/turn  = ",
     +   2.00D0*pi/avtheta," +/- ",2.0D0*pi*sdtheta/(avtheta**2)
c
        enddo 
       endif 
c
       call exit(0)
       stop
       end


       subroutine arcdist(a,b,c,d,i,j,k,l,m,n,R)
c
c      Determine the distance of a point (i,j,k) from the 
c      line of intersection of a sphere (with centre at (a,b,c)
c      and radius R ie. (x-a)**2 + (y-b)**2 + (z-c)**2 = R)
c      and a plane (lying on the point (a,b,c) 
c      ie. l*(x-a) + m*(y-b) + n*(z-c) = 0, where {l,m,n,R}>0).
c
       real*8    a,b,c,i,j,k,l,m,n,R
       real*8    d,d2,tmp1,tmp2,tmp3,tmp4
c
       tmp1=((i-a)**2)+((j-b)**2)+((k-c)**2)
       tmp2=((l*(i-a))+(m*(j-b))+(n*(k-c)))**2
       tmp3=(l*l)+(m*m)+(n*n)
       tmp4=tmp2/tmp3
       d=0.0D0
       d2=tmp1+(R*R)-(2.0D0*R*dsqrt(tmp1-tmp4))
       d=dsqrt(d2)
c
       return
       end


       subroutine arc(stage)
c
       real*8    x(500),y(500),z(500)
       real*8    u(500),v(500),w(500)
       real*8    a,b,c,d,l,m,n,R
       real*8    dd2da,dd2db,dd2dc,dd2dl,dd2dm,dd2dn,dd2dR
       real*8    ds2da,ds2db,ds2dc,ds2dl,ds2dm,ds2dn,ds2dR
       real*8    param(7),grad(7)
       real*8    dav,d2,s2
       real*8    objective
       real*8    tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7
       integer*4 natoms
       integer*4 stage
       integer*4 fix_radius
       integer*4 i
c
       common/coord/x,y,z,u,v,w,natoms
       common/cg/objective,param,grad,fix_radius
c
c      Substitute parameter names
c
       a=param(1)
       b=param(2)
       c=param(3)
       l=param(4)
       m=param(5)
       n=param(6)
       R=param(7)
c
c      Initialize variables
c
       dav=0.0D0
       d2=0.0D0
       s2=0.0D0
       dd2da=0.0D0
       dd2db=0.0D0
       dd2dc=0.0D0
       dd2dl=0.0D0
       dd2dm=0.0D0
       dd2dn=0.0D0
       dd2dR=0.0D0
       ds2da=0.0D0
       ds2db=0.0D0
       ds2dc=0.0D0
       ds2dl=0.0D0
       ds2dm=0.0D0
       ds2dn=0.0D0
       ds2dR=0.0D0
c
c      Determine Sum(d**2) and the variance - the function to minimize
c
       do i=1,natoms
        call arcdist(a,b,c,d,x(i),y(i),z(i),l,m,n,R)
        dav=dav+d
        d2=d2+d*d
       enddo
       dav=dav/float(natoms)
       do i=1,natoms
        call arcdist(a,b,c,d,x(i),y(i),z(i),l,m,n,R)
        s2=s2+(d-dav)*(d-dav)
       enddo
       s2=s2/float(natoms-1)
c
       if(stage.eq.0) then
c
c      Determine derivatives of Sum(d**2) wrt each variable (a,b,c,l,m,n,R)
c
       tmp1=(l*l)+(m*m)+(n*n)
       do i=1,natoms
        tmp2=(l*(x(i)-a))+(m*(y(i)-b))+(n*(z(i)-c))
        tmp3=tmp2/tmp1
        tmp4=tmp2*tmp3
        tmp5=tmp3*tmp3
        tmp6=((x(i)-a)**2)+((y(i)-b)**2)+((z(i)-c)**2)
        tmp7=tmp6-tmp4
c
        dd2da=dd2da+2.0D0*(a-x(i))-2.0D0*R*(a-x(i)+l*tmp3)/dsqrt(tmp7)
        dd2db=dd2db+2.0D0*(b-y(i))-2.0D0*R*(b-y(i)+m*tmp3)/dsqrt(tmp7)
        dd2dc=dd2dc+2.0D0*(c-z(i))-2.0D0*R*(c-z(i)+n*tmp3)/dsqrt(tmp7)
        dd2dl=dd2dl+2.0D0*R*(((x(i)-a)*tmp3)-(l*tmp5))/dsqrt(tmp7)
        dd2dm=dd2dm+2.0D0*R*(((y(i)-b)*tmp3)-(m*tmp5))/dsqrt(tmp7)
        dd2dn=dd2dn+2.0D0*R*(((z(i)-c)*tmp3)-(n*tmp5))/dsqrt(tmp7)
        dd2dR=dd2dR+(2.0D0*R)-(2.0D0*dsqrt(tmp7))
       enddo
c
c      Substitute the gradient array
c
       grad(1)=dd2da
       grad(2)=dd2db
       grad(3)=dd2dc
       grad(4)=dd2dl
       grad(5)=dd2dm
       grad(6)=dd2dn
       grad(7)=dd2dR
       objective=d2
c
       else
c
c      Determine derivatives of variance (s**2) wrt each variable (a,b,c,l,m,n,R)
c
       tmp1=(l*l)+(m*m)+(n*n)
       do i=1,natoms
        tmp2=(l*(x(i)-a))+(m*(y(i)-b))+(n*(z(i)-c))
        tmp3=tmp2/tmp1
        tmp4=tmp2*tmp3
        tmp5=tmp3*tmp3
        tmp6=((x(i)-a)**2)+((y(i)-b)**2)+((z(i)-c)**2)
        tmp7=tmp6-tmp4
        call arcdist(a,b,c,d,x(i),y(i),z(i),l,m,n,R)
        dd2da=2.0D0*(a-x(i))-2.0D0*R*(a-x(i)+l*tmp3)/dsqrt(tmp7)
        ds2da=ds2da+(1.0D0-(dav/d))*dd2da
        dd2db=2.0D0*(b-y(i))-2.0D0*R*(b-y(i)+m*tmp3)/dsqrt(tmp7)
        ds2db=ds2db+(1.0D0-(dav/d))*dd2db
        dd2dc=2.0D0*(c-z(i))-2.0D0*R*(c-z(i)+n*tmp3)/dsqrt(tmp7)
        ds2dc=ds2dc+(1.0D0-(dav/d))*dd2dc
        dd2dl=2.0D0*R*(((x(i)-a)*tmp3)-(l*tmp5))/dsqrt(tmp7)
        ds2dl=ds2dl+(1.0D0-(dav/d))*dd2dl
        dd2dm=2.0D0*R*(((y(i)-b)*tmp3)-(m*tmp5))/dsqrt(tmp7)
        ds2dm=ds2dm+(1.0D0-(dav/d))*dd2dm
        dd2dn=2.0D0*R*(((z(i)-c)*tmp3)-(n*tmp5))/dsqrt(tmp7)
        ds2dn=ds2dn+(1.0D0-(dav/d))*dd2dn
        dd2dR=(2.0D0*R)-(2.0D0*dsqrt(tmp7))
        ds2dR=ds2dR+(1.0D0-(dav/d))*dd2dR
       enddo
c
c      Substitute the gradient array
c
       grad(1)=ds2da/float(natoms-1)
       grad(2)=ds2db/float(natoms-1)
       grad(3)=ds2dc/float(natoms-1)
       grad(4)=ds2dl/float(natoms-1)
       grad(5)=ds2dm/float(natoms-1)
       grad(6)=ds2dn/float(natoms-1)
       grad(7)=ds2dR/float(natoms-1)
       objective=s2
c
       endif
c
       return
       end


       subroutine acgm(stage,tol,ier,iter,step)
c
c     A Conjugate Gradients Minimizer
c
c     The procedure adopted here follows that described in the following...
c
c     "Function minimization by conjugate gradients"
c      R. Fletcher and C.M. Reeves.
c      British Computer Journal Vol. 7, pp. 149-154 (1964)
c
c     "Restart procedures for the conjugate gradient method"
c      M.J.D. Powell
c      Mathematical Programming Vol. 12, pp. 241-254 (1977)
c
c     "Conjugate-gradient methods for large-scale minimization in meteorology"
c      I.M. Navon and D.M. Legler
c      Monthly Weather Review Vol. 115, pp. 1479-1501 (1987)
c
       real*8     objective,param(7),grad(7)
       real*8     param_min(7),grad_min(7),search(7)
       real*8     search_rstrt(7),grad_rstrt(7),g0(7),gg(7)
       real*8     a,b,te,beta,gamma,step
       real*8     dotpga,dotpgb,dotgaga,dotgbgb,dotpg0,dotgbg0
       real*8     dotdgg,dotg0g0,dotgbga
       real*8     fmin,fold
       real*8     small,gmax
       real*8     tol
       integer*4  stage,ier,iter
       integer*4  i,j
       integer*4  maxit,it,xit,mit,kmax,kount,ik,ns
       integer*4  Halve_step,Halve_step_max
       integer*4  fix_radius
       integer*4  out
c
       common/cg/objective,param,grad,fix_radius
       common/io/out
c
       maxit=100000  !  Absurdly large number of iterations
       mit=100  !  Absurdly large number of attempts at quadratic interpolation
       small=1.0D-64
       kmax=16
       Halve_step_max=64
c
       call arc(stage)
       fmin=objective
       dotpga=0.0D0
       dotgaga=0.0D0
       do i=1,(6+fix_radius)
        search(i)=(-grad(i))
        grad_rstrt(i)=grad(i)
        search_rstrt(i)=search(i)
        param_min(i)=param(i)
        grad_min(i)=grad(i)
        dotpga=dotpga+(search(i)*grad(i))
        dotgaga=dotgaga+(grad(i)*grad(i))
        g0(i)=grad(i) 
       enddo
       a=0.0D0
       gmax=-99999.99999D0
       do i=1,(6+fix_radius)
        gmax=max(gmax,grad(i))
       enddo
       if(stage.eq.0) then
        step=gmax*gmax/dotgaga
       endif
       b=step
       beta=0.0D0
       gamma=0.0D0
       iter=0
       ik=1
       it=1
       ns=1
c
c     Start the iterations
c
       do j=1,maxit
c
        dotpg0=0.0D0
c
c     Ensure the radius stays positive
c
        Halve_step=0
        do while (((param_min(7)+(step*search(7))).le.0.0D0).and.
     +             (Halve_step.lt.Halve_step_max))
         step=step/2.0D0
         Halve_step=Halve_step+1
        enddo
c
        do i=1,(6+fix_radius)
         param(i)=param_min(i)+(step*search(i))
        enddo
        call arc(stage)
        dotpgb=0.0D0
        do i=1,(6+fix_radius)
         dotpgb=dotpgb+(grad(i)*search(i))
        enddo
c
c     Interpolate if the objective increases (else, extrapolate).
c
        if((objective.ge.fmin)) then
c
        xit=0
 1      continue
c
c     Quadratic interpolation.
c     Note, te is the estimate of the minimum lying between a and b.
c
       if((abs(a-b).lt.small).or.(abs(dotpga-dotpgb).lt.small)) then
        fold=fmin
        goto 2
       endif
c
       te=((a**2-b**2)/(2.0D0*(a-b)))-((fmin-objective)/(dotpga-dotpgb))
c
        do i=1,(6+fix_radius)
         param(i)=param_min(i)+(te*search(i))
        enddo
        fold=objective
        call arc(stage)
c
c     Check that interpolation led to a decrease in the objective,
c     or that the decrease in the objective is sufficiently small.
c
 2      continue
        if(((objective.lt.fold).and.(objective.lt.fmin)).or.
     +      (abs(fold-fmin).lt.tol)) then
c
         dotgaga=0.0D0
         dotgbga=0.0D0
         dotgbgb=0.0D0
         dotgbg0=0.0D0
         dotpga=0.0D0
         dotpgb=0.0D0
         dotpg0=0.0D0
         dotg0g0=0.0D0
         do i=1,(6+fix_radius)
          dotgaga=dotgaga+(grad_min(i)*grad_min(i))
          dotgbga=dotgbga+(grad(i)*grad_min(i))
          dotgbgb=dotgbgb+(grad(i)*grad(i))
          dotgbg0=dotgbg0+(grad(i)*g0(i))
          dotpga=dotpga+(search(i)*grad_min(i))
          dotpgb=dotpgb+(search(i)*grad(i))
          dotpg0=dotpg0+(search(i)*g0(i))
          dotg0g0=dotg0g0+(g0(i)*g0(i))
          param_min(i)=param(i)
          grad_min(i)=grad(i)
         enddo
c
c     Test for convergence.
c
         dotgaga=0.0D0
         do i=1,(6+fix_radius)
          dotgaga=dotgaga+(grad(i)*grad(i))
         enddo
         if(dotgaga.le.tol) then
          ier=0
          return
         endif
c
c     Update beta (a choice of methods - which to choose?)
c
         if(abs(dotpgb-dotpg0).gt.small) then 
c         beta=(dotgbgb)/(dotg0g0)  !  Fletcher-Reeves
c         beta=(dotgbgb-dotgbg0)/(dotg0g0)  !  Polak-Ribiere
          beta=(dotgbgb-dotgbg0)/(dotpgb-dotpg0)  !  Shanno-Perry
         endif
c
c     Update Beale's gamma (if required)
c
         if((abs(dotdgg).gt.small).and.(ik.gt.it+1)) then
          gamma=0.0D0
          do i=1,(6+fix_radius)
           gamma=gamma+gg(i)*grad(i)
          enddo
          gamma=gamma/dotdgg
         endif 
c
         fmin=objective
         do i=1,(6+fix_radius)
          g0(i)=grad(i)
         enddo
         ik=ik+1
         ns=ns+1
c
c     After N iterations, where N equals the number of parameters, 
c     either reset the search direction to steepest descents, 
c      .OR.
c     set t=k-1 to signify a restart is warranted.
c
         if(ns.ge.(6+fix_radius)) then
c         gamma=0.0D0
c         beta=0.0D0
          it=ik-1
          ns=0
c
c     Test Powell's inequality (to determine if a restart is necessary)
c
         elseif(abs(dotgbg0).ge.0.2D0*dotgbgb) then
c
c     When there is very little orthogonality remaining between the 
c     current and previous gradients, either the search direction is 
c     reset to steepest descents,
c      .OR.
c     set t=k-1 to signify a restart is warranted.
c
c         gamma=0.0D0
c         beta=0.0D0
          it=ik-1
c
c     Restart the search direction if the search direction is not
c     sufficiently downhill.
c
         elseif((dotpgb.gt.(-0.8D0)*dotgbgb).or.
     +          (dotpgb.lt.(-1.2D0)*dotgbgb)) then
          it=ik-1
         endif
c
c     Reset the restart direction if t=k-1, and determine 
c     Beale's parameters
c
         if(it.eq.ik-1) then
          do i=1,(6+fix_radius)
           search_rstrt(i)=search(i) 
          enddo 
          dotdgg=0.0D0 
          do i=1,(6+fix_radius)
           gg(i)=grad(i)-grad_rstrt(i)
           grad_rstrt(i)=grad(i)
           dotdgg=dotdgg+(search_rstrt(i)*gg(i))
          enddo 
         endif
c
c     New search direction
c
         dotpga=0.0D0
         do i=1,(6+fix_radius)
c
c     Beale's gamma set to zero if k = t+1 (recall, k > t)
c
          if(it.eq.ik-1) then
           search(i)=(-grad(i))+(beta*search(i))
          else
           search(i)=(-grad(i))+(beta*search(i))+(gamma*search_rstrt(i))
          endif
          dotpga=dotpga+(search(i)*grad(i))
         enddo
c
c     Force the search direction to be downhill.
c
         if(dotpga.gt.0.0D0) then
          dotpga=0.0D0
          do i=1,(6+fix_radius)
           search(i)=(-grad(i))
           dotpga=dotpga+(search(i)*grad(i))
          enddo
         endif
c
c     Determine the new step length.
c
         gmax=-99999.99999
         do i=1,(6+fix_radius)
          gmax=max(gmax,grad(i))
         enddo
         step=gmax*gmax/dotgaga
         a=0.0D0
         b=step
c
        else
c
c     Interpolation did not locate a minimum. Continue the search over
c     the intervals a-te, te-b, or back-track (b to a) in half steps.
c
         if(xit.gt.mit) then
c
c         "Too many attempts at quadratic interpolation"
c
          ier=1
          return
         endif
c
         if((objective.lt.fold).and.(objective.gt.fmin)) then
c
c     a-te interpolation.
c
          b=te
          dotpgb=0.0D0
          do i=1,(6+fix_radius)
           dotpgb=dotpgb+(search(i)*grad(i))
          enddo
          xit=xit+1
          goto 1
         elseif((objective.gt.fold).and.(objective.lt.fmin)) then
c
c     te-b interpolation.
c
          a=te
          dotpga=0.0D0
          do i=1,(6+fix_radius)
           dotpga=dotpga+(search(i)*grad(i))
          enddo
          fmin=objective
          objective=fold
          xit=xit+1
          goto 1
         else
c
c     Reduce the step by halves until the objective is less than
c     the previous best.
c
          kount=0
          do while((objective.gt.fmin).and.(kount.lt.kmax))
           b=(a+b)/2.0D0
           do i=1,(6+fix_radius)
            param(i)=param_min(i)+(b*search(i))
           enddo
           call arc(stage)
           kount=kount+1
          enddo
          if(kount.lt.kmax) then
           fmin=objective
           do i=1,(6+fix_radius)
            param_min(i)=param(i)
            grad_min(i)=grad(i)
           enddo
          endif
          fold=fmin
          xit=xit+1
          goto 2
c
         endif
        endif
c
        else
c
c     Force a restart if changes in step or gradient are too small.
c
       if((abs(a-b).lt.small).or.(abs(dotpga-dotpgb).lt.small)) then
        fold=fmin
        goto 2
       endif
c
c     Quadratic extrapolation
c
       te=((a**2-b**2)/(2.0D0*(a-b)))-((fmin-objective)/(dotpga-dotpgb))
c
c     Force a restart if the search direction changes.
c
        if((te*(b-a)).lt.0.0D0) then
         fold=fmin
         goto 2
        endif
c
         a=0.0D0
         fmin=objective
         dotpga=0.0D0
         dotgaga=0.0D0
         do i=1,(6+fix_radius)
          param_min(i)=param(i)
          grad_min(i)=grad(i)
          dotpga=dotpga+(search(i)*grad(i))
         enddo
         b=a+te
         step=te
        endif
c
        iter=iter+1
       enddo   !  End of main loop
c
       if(iter.ge.maxit) then
c
c       "Too many steps in conjugate gradients minimization"
c
        ier=2
        return
       endif
c
       return
       end


       subroutine UVW(a,b,c,i,j,k,l,m,n,R,u,v,w)
c
       real*8    a,b,c,i,j,k,l,m,n,R
       real*8    up,vp,wp,u,v,w
       real*8    tmp1,tmp2,tmp3
c
       tmp1=(l*a)+(m*b)+(n*c)
       tmp2=((l*i)+(m*j)+(n*k)-tmp1)/((l*l)+(m*m)+(n*n))
       up=i-l*tmp2
       vp=j-m*tmp2
       wp=k-n*tmp2
       tmp3=dsqrt(((a-up)*(a-up))+((b-vp)*(b-vp))+((c-wp)*(c-wp)))
       u=a+R*(up-a)/tmp3
       v=b+R*(vp-b)/tmp3
       w=c+R*(wp-c)/tmp3
c
       return
       end


       subroutine res2res(resid,numer,char)
c
c      Unpack the residue identifier
c
       integer*4    numer,i,j
       character*80 resid
       character*1  char
c
       char(1:1)=' '
       j=80
c
       do i=80,1,-1
        if((ichar(resid(i:i)).ge.ichar('A')).and.
     +     (ichar(resid(i:i)).le.ichar('Z'))) then
         char(1:1)=resid(i:i)
         j=i-1
        elseif(resid(i:i).eq.' ') then
         read(resid(i+1:j),*)numer
         goto 1
        endif
       enddo
c
 1     continue
       return
       end


       subroutine rhj(String)
c
c      Right-hand justify a string (of 80 characters)
c
       character*80 String
       integer*4 i,j
c
       do i=1,80
        if(String(80:80).eq.' ') then
         do j=79,1,-1
          String(j+1:j+1)=String(j:j)
         enddo
         String(1:1)=' '
        endif
       enddo
c
       return
       end


       subroutine lhj(String)
c
c      Left-hand justify a string (of 4 characters)
c
       character*4 String
       integer*4 i,j
c
       do i=1,4
        if(String(1:1).eq.' ') then
         do j=1,3
          String(j:j)=String(j+1:j+1)
         enddo
         String(4:4)=' '
        endif
       enddo
c
       return
       end

       subroutine initialize
c 
       real*8     x(500),y(500),z(500)
       real*8     u(500),v(500),w(500)
       real*8     objective,param(7),grad(7)
       real*8     locus(3),mid(3),a(3),q(3),e(3)
       real*8     small,d,offset,norm,t
       integer*4  natoms
       integer*4  fix_radius
       integer*4  i
c
       common/coord/x,y,z,u,v,w,natoms
       common/cg/objective,param,grad,fix_radius
       common/io/out
c
       small=1.0D-64
c
c      Identify atom locus
c
       locus(1)=0.0D0
       locus(2)=0.0D0
       locus(3)=0.0D0
       do i=1,natoms
        locus(1)=locus(1)+x(i)
        locus(2)=locus(2)+y(i)
        locus(3)=locus(3)+z(i)
       enddo
       do i=1,3
        locus(i)=locus(i)/float(natoms)
       enddo
c
c      Identify mid-point on the line between first and last points
c
       mid(1)=(x(natoms)+x(1))/2.0D0
       mid(2)=(y(natoms)+y(1))/2.0D0
       mid(3)=(z(natoms)+z(1))/2.0D0
c
c      Place great circle origin a distance 'offset' from mid-point 
c
       a(1)=x(natoms)-x(1)
       a(2)=y(natoms)-y(1)
       a(3)=z(natoms)-z(1)
       norm=sqrt(a(1)**2 + a(2)**2 + a(3)**2)
       if(norm.le.small) then
        write(6,'(A)')"Error initializing CG parameters"
        write(6,'(A)')"Program must stop"
        stop
       endif
       do i= 1,3
        a(i)=a(i)/norm
       enddo
       e(1)=((((locus(1)-x(1))*a(1))+
     +       ((locus(2)-y(1))*a(2))+
     +       ((locus(3)-z(1))*a(3)))*a(1)) -
     +       (locus(1)-x(1))
       e(2)=((((locus(1)-x(1))*a(1))+
     +       ((locus(2)-y(1))*a(2))+
     +       ((locus(3)-z(1))*a(3)))*a(2)) -
     +       (locus(2)-y(1))
       e(3)=((((locus(1)-x(1))*a(1))+
     +       ((locus(2)-y(1))*a(2))+
     +       ((locus(3)-z(1))*a(3)))*a(3)) -
     +       (locus(3)-z(1))
       norm=sqrt(e(1)**2 + e(2)**2 + e(3)**2)
       if(norm.le.small) then
        do i=1,3
         if(a(i).ne.0.0D0) then
          e(i)=-1.0D0/a(i)
         else
          e(i)=0.0D0
         endif
        enddo
       else
        do i=1,3
         e(i)=e(i)/norm
        enddo
       endif
c
c      Initial radius of curvature is either
c       1. provided as a seed,
c      or
c       2. estimated, with a minimum of 20.0, a maximum of 1000.0,
c          and scaling with the displacement of the locus from 
c          the line between first and last points.
c
       if(param(7).eq.0.0D0) then
        param(7)=20.0D0 + max(((10.0D0-norm)/10.0D0),0.0D0)*980.0D0
       endif
c
c      Determine offset of great circle origin from mid-point
c
       d=sqrt((mid(1)-x(1))**2 +
     +        (mid(2)-y(1))**2 +
     +        (mid(3)-z(1))**2)
       if(d.lt.param(7)) then
        offset=sqrt(param(7)**2-d**2)
       else
        offset=param(7)
       endif
       do i=1,3
        param(i)=mid(i) + offset*e(i)
       enddo
c
c      Plane orientation MUST be non-zero
c
       param(4)=y(1)*z(natoms) - param(3)*z(natoms) -
     +          y(1)*param(3) - y(natoms)*z(1) +
     +          y(natoms)*param(3) + param(2)*z(1)
       param(5)=x(1)*z(natoms) - param(1)*z(natoms) -
     +          x(1)*param(3) - x(natoms)*z(1) +
     +          x(natoms)*param(3) + param(1)*z(1)
       param(6)=x(1)*y(natoms) - param(1)*y(natoms) -
     +          x(1)*param(2) - x(natoms)*y(1) +
     +          x(natoms)*param(2) + param(1)*y(1)
c
       norm=sqrt(param(4)**2 +param(5)**2 + param(6)**2)
       if(norm.lt.small) then
        write(6,'(A)')"Error initializing CG parameters"
        write(6,'(A)')"Program must stop"
        stop
       else
        do i=4,6
         param(i)=param(i)/norm
        enddo
       endif
c
       return
       end

       subroutine helixdir(dir)
c
       real*8    x(500),y(500),z(500)
       real*8    u(500),v(500),w(500)
       real*8    a,b,c,d,l,m,n,R
       real*8    param(7),grad(7)
       real*8    objective
       real*8    a1,a2,a3
       real*8    norm
       real*8    c1,c2,c3
       real*8    cdota
       real*8    cparalell1,cparalell2,cparalell3
       real*8    cperp1,cperp2,cperp3
       real*8    dir(6)
       integer*4 natoms
       integer*4 fix_radius
       integer*4 i,j
c
       common/coord/x,y,z,u,v,w,natoms
       common/cg/objective,param,grad,fix_radius
c
       j=1
       c1=u(natoms)-u(1)
       c2=v(natoms)-v(2)
       c3=w(natoms)-w(3)
c
       do i=1,natoms,natoms-1
        a1=u(i)-param(1)
        a2=v(i)-param(2)
        a3=w(i)-param(3)
        norm=dsqrt(a1**2 + a2**2 + a3**2)
        a1=a1/norm
        a2=a2/norm
        a3=a3/norm
c 
        cdota=c1*a1 + c2*a2 + c3*a3
        cparallel1=cdota*a1
        cparallel2=cdota*a2
        cparallel3=cdota*a3
c
        cperp1=c1-cparallel1
        cperp2=c2-cparallel2
        cperp3=c3-cparallel3
        norm=dsqrt(cperp1**2 + cperp2**2 + cperp3**2)
        dir(j)=cperp1/norm
        dir(j+1)=cperp2/norm
        dir(j+2)=cperp3/norm
        j=j+3
       enddo
c
       return
       end