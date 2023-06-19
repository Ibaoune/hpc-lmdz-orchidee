subroutine dates_demo
! --------------------------------------------------------------
!
! Conseils a l'utilisateur:
!
! 1. VOUS COMPILEZ LES ENTIERS EN 32 BITS:
! Utilisez alors les routines
! - ecartds: Ecart en secondes entre deux dates.
! - ecartdj: Ecart en jours entre deux dates.
! - dapluss: Date dans n secondes.
! - daplusj: Date dans n jours.
! - qqmmaa: Conversion d'un entier type AAAAQQMM vers une date en clair.
! - ijoursem: Jour de la semaine de la date d'entree.
! - quant: quantieme de l'annee d'une date donnee.
! Ces routines sont compatibles avec des entiers 32 bits.
! En effet elles appelent les routines citees ci-dessous, mais avec
! les parametres subsequents assurant que seuls des entiers
! representables en 32 bits y soient utilises.
!
! 2. VOUS COMPILEZ LES ENTIERS EN 64 BITS:
! Vous pouvez alors utiliser toutes les routines ci-dessus
! plus les suivantes, qui traitent des formats de dates
! en entree/sortie en JOURS, HEURES, MINUTES ou SECONDES:
! - ecartd: Ecart entre deux dates.
! - gregod: Conversion Date > Ecart par rapport a une date fixe.
! - gregoi: Conversion Ecart par rapport a une date fixe > Date.
! - daplus: Quelle sera la date dans n jours (ou heures, etc...)?
! - amqhms_vers_dj: Conversion date grégorienne (en 5 entiers et un réel) > date julienne.
! - dj_vers_amqhms: Conversion date julienne > date grégorienne (en 5 entiers et un réel).
! - amqhmsree_vers_dj: Conversion date grégorienne (en un seul réel) > date julienne.
! - dj_vers_amqhmsree: Conversion date julienne > date grégorienne (en un seul réel).
!
! --------------------------------------------------------------
!
! Définition des dates employées ci-dessous:
!
! Date julienne DJ:
!       Elle est composée d'un réel.
!       R1: Ce réel croît de 1 tous les jours,
!               et vaut 2451545.0 le 1er janvier 2000 à 12 UTC.
!
! Date grégorienne "en clair" AMQHMS:
!       Elle est composée de 5 entiers et d'un réel.
!       E1: Année (4 chiffres!)
!       E2: Mois
!       E3: Jour
!       E4: Heure
!       E5: Minute
!       R1: Seconde
! --------------------------------------------------------------


IMPLICIT NONE
end
subroutine date_plus_ech(kan,kmo,kqu,psssss,pstati,cdtit)
! --------------------------------------------------------------
! Ecriture en clair d'une date de type BASE 2000.01.15 00:00 +72H VALID 2000.01.18 15:00.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! kan,kmo,kqu,psssss,pstati
! En sortie:
! cdtit
! --------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: kan,kmo,kqu,ihe,imi,imiv,ihev,iquv,imov,ianv,ilze
REAL(KIND=8) :: psssss,pstati
REAL(KIND=8) :: zs
REAL(KIND=8) :: zsssss,zdj,zsv
REAL(KIND=8) :: zech
character*200 clzue,clze,clech
character *(*) cdtit
!
!-------------------------------------------------
! Date de validité.
!-------------------------------------------------
!
zs=0.
zsssss=psssss/3600.
ihe=int(zsssss) ! heure de la base.
zsssss=(zsssss-real(ihe))*60.
imi=int(zsssss) ! minute de la base.
zsssss=zsssss-real(imi)
call amqhms_vers_dj(kan,kmo,kqu,ihe,imi,zs,zdj)
zdj=zdj+pstati/86400. ! date julienne de validité.
call dj_vers_amqhms(zdj,ianv,imov,iquv,ihev,imiv,zsv) ! date grégorienne de validité.
if(pstati < 3600.) then
!
!-------------------------------------------------
! Echéance en minutes.
!-------------------------------------------------
!
	zech=pstati/60. ; clzue='mn'
elseif(pstati < 259200.) then
!
!-------------------------------------------------
! Echéance en heures.
!-------------------------------------------------
!
	zech=pstati/3600. ; clzue='h'
else
!
!-------------------------------------------------
! Echéance en jours.
!-------------------------------------------------
!
	zech=pstati/86400. ; clzue='j'
endif
!
! Affichage de l'echeance avec deux chiffres apres la virgule.
!
write(clze,fmt='(f9.2)') zech
!
! Si l'echeance est voisine d'un entier a mieux que 10**-2 pres,
! on l'affiche au format entier.
!
if(clze(len_trim(clze)-2:len_trim(clze)) == '.00') then
	clze=clze(1:len_trim(clze)-3)
endif
clze=adjustl(clze)
ilze=len_trim(clze)
clech=clze(1:ilze)//clzue
!
!-------------------------------------------------
! Titre 3, de type
! BASE 2000.01.15 00:00 +72H VALID 2000.01.18 15:00.
!-------------------------------------------------
!
if(imi == 0 .and. imiv == 0) then
!
!-------------------------------------------------
! Les minutes de base et validité sont nulles.
! On ne les affiche pas.
!-------------------------------------------------
!
	write(cdtit,fmt='(a,i2,a,i2.2,a,i4.4,a,i2.2,3a,i2,a,i2.2,a,i4.4,a,i2.2,a)')&
	&'BASE ',kqu,'.',kmo,'.',kan,' ',ihe,'h UTC + ',clech(1:len_trim(clech))&
	&,', VALID ',iquv,'.',imov,'.',ianv,' ',ihev,'h UTC'
else
!
!-------------------------------------------------
! Les minutes de base ou validité sont non nulles.
! On les affiche.
!-------------------------------------------------
!
	write(cdtit,fmt='(a,i2,a,i2.2,a,i4.4,a,i2.2,a,i2.2,3a,i2,a,i2.2,a,i4.4,a,i2.2,a,i2.2,a)')&
	&'BASE ',kqu,'.',kmo,'.',kan,' ',ihe,':',imi,' UTC + ',clech(1:len_trim(clech))&
	&,' VALID ',iquv,'.',imov,'.',ianv,' ',ihev,':',imiv,' UTC'
endif
end
subroutine datc(kaaaa,kmm,kqq,khh,kmi,kss,kjs,cdjs,cddt)
! --------------------------------------------------------------
! **** *datc* Date courante machine.
! --------------------------------------------------------------
! Sujet:
! ------
! Arguments explicites:
! ---------------------
! Arguments implicites:
! ---------------------
! Methode:
! --------
! Externes:
! ---------
! Auteur:   95-05, J.M. Piriou.
! -------
! Modifications:
! --------------------------------------------------------------
! En entree:
! En sortie:
! kaaaa      annee.
! kmm      mois.
! kqq      quantieme.
! khh      heure.
! kmi      minute.
! kss      seconde.
! kjs      jour de la semaine (0: dimanche, 6 samedi).
! cdjs      jour de la semaine sur 3 caracteres (Dim, Lun, etc...).
! cddt      date totale (19950301-Mer-16:56:32).
! --------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: idatat(8)
INTEGER(KIND=4) :: kjs
INTEGER(KIND=4) :: kss
INTEGER(KIND=4) :: kmi
INTEGER(KIND=4) :: khh
INTEGER(KIND=4) :: kqq
INTEGER(KIND=4) :: kmm
INTEGER(KIND=4) :: kaaaa
INTEGER(KIND=4) :: iaaaammqq
INTEGER(KIND=4) :: ijoursem
REAL(KIND=8) :: zs
character*200 clgol1,clgol2,clgol3
character*3 cdjs
character*(*) cddt
character*3 cljour(0:6)
data cljour/'Dim','Lun','Mar','Mer','Jeu','Ven','Sam'/
!
!-------------------------------------------------
! Date courante à la f90.
!-------------------------------------------------
!
clgol1=' '
clgol2=' '
clgol3=' '
call date_and_time(clgol1,clgol2,clgol3,idatat)
!
!-------------------------------------------------
! clgol1 est du type "AAAAMMQQ".
!-------------------------------------------------
!
read(clgol1,fmt='(i4,2i2)') kaaaa,kmm,kqq
!
!-------------------------------------------------
! clgol2 est du type "HHMMSS.SSS".
!-------------------------------------------------
!
read(clgol2,fmt='(2i2)') khh,kmi
read(clgol2(5:),fmt=*) zs
kss=nint(zs)
read(clgol1,fmt='(i8)') iaaaammqq
!
!-------------------------------------------------
! Jour de la semaine.
!-------------------------------------------------
!
kjs=ijoursem(iaaaammqq)
cdjs=cljour(kjs)
!
!-------------------------------------------------
! Date totale.
!-------------------------------------------------
!
write(cddt,fmt='(i4.4,a,2(i2.2,a),2a,i2.2,a,i2.2,a,i2.2)') &
&kaaaa,'_',kmm,'_',kqq,'_',cdjs,'_',khh,':',kmi,':',kss
end
subroutine amqhms_vers_dj(kaaaa,kmm,kqq,khh,kmn,ps,pdj)
! --------------------------------------------------------------------------
! **** *amqhms_vers_dj*
! --------------------------------------------------------------------------
! Auteur:
! -------
! 1999-08-17, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! kaaaa année (4 chiffres!)
! kmm   mois
! kqq   quantième du mois
! khh   heure
! kmn   minute
! ps    seconde
! En sortie:
! pdj date julienne associée à la date grégorienne UTC d'entrée
! --------------------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: IDATE1
INTEGER(KIND=4) :: IDATE2
INTEGER(KIND=4) :: IECART
INTEGER(KIND=4) :: KAAAA
INTEGER(KIND=4) :: KHH
INTEGER(KIND=4) :: KMM
INTEGER(KIND=4) :: KMN
INTEGER(KIND=4) :: KQQ
REAL(KIND=8) :: PDJ
REAL(KIND=8) :: PS

idate1=20000101
idate2=kaaaa*10000+kmm*100+kqq
!
!-------------------------------------------------
! Nombre de jours écoulés entre la date
! d'entrée à 0h UTC et le 1er janvier 2000 à 0h UTC.
!-------------------------------------------------
!
call ecartdj(idate1,idate2,iecart)
!
!-------------------------------------------------
! Date julienne.
!-------------------------------------------------
!
pdj=2451545.0- 0.5 +real(iecart)+real(khh)/24. &
& +real(kmn)/1440.+ps/86400.
end
subroutine daplus(kdat1,kopt,kdelt,kdat2)
! --------------------------------------------------------------------------
! **** *DAPLUS* Quelle sera la date dans n jours (ou heures, etc...)?
! --------------------------------------------------------------------------
! Auteur:
! -------
! 94-10-31, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! kdat1
! kopt option de precision sur les dates:
! 1 : au jour pres
! 2 : a l'heure pres
! 3 : a la minute pres
! 4 : a la seconde pres
! - si kopt=1 : kdat au format AAAAMMQQ
! - si kopt=2 : kdat au format AAAAMMQQHH
! - si kopt=3 : kdat au format AAAAMMQQHHMM
! - si kopt=4 : kdat au format AAAAMMQQHHMMSS
! (cf. GREGOD).
! kdelt duree a ajouter a kdat1, unite: celle imposee par kopt.
! En sortie:
! kdat2 date finale.
!
! --------------------------------------------------------------------------
! Exemple: call DAPLUS(19940503,1,456,ires) fournira
! dans ires la date au format AAAAMMQQ situee 456 jours apres
! le 3 mai 1994.
! --------------------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! 1. Vous compilez les entiers sur 32 bits:
! Vous devez alors vous limiter a kopt <= 2.
! 2. Vous compilez les entiers sur 64 bits:
! Vous pouvez utiliser toutes les valeurs de kopt.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: IGRE
INTEGER(KIND=4) :: KDAT1
INTEGER(KIND=4) :: KDAT2
INTEGER(KIND=4) :: KDELT
INTEGER(KIND=4) :: KOPT
call gregod(kdat1,kopt,igre)
igre=igre+kdelt
call gregoi(igre,kopt,kdat2)
end
subroutine daplusj(k1,kec,k2)
! --------------------------------------------------------------
! **** *daplusj* Date dans n jours.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   97-11, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! k1 date de depart au format AAAAMMQQ.
! kec nombre de jours ecoules.
! En sortie:
! k2 date d'arrivee au format AAAAMMQQ.
! --------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! PRECISION:
! Cette routine est utilisable avec des entiers 32 bits ou 64 bits.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! -------------------------------------------------
! Date d'arrivee au jour pres.
! -------------------------------------------------
!


IMPLICIT NONE
INTEGER(KIND=4) :: K1
INTEGER(KIND=4) :: K2
INTEGER(KIND=4) :: KEC
call daplus(k1,1,kec,k2)
end
subroutine dapluss(cd1,kec,cd2)
! --------------------------------------------------------------
! **** *dapluss* Date dans n secondes.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   97-11, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! cd1 date de depart au format 'AAAAMMQQHHNNSS'.
! kec nombre de secondes ecoulees.
! En sortie:
! cd2 date d'arrivee au format 'AAAAMMQQHHNNSS'.
! --------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! Cette routine est utilisable avec des entiers 32 bits,
! si l'ecart entre les deux dates est inferieur a 2**31 secondes,
! soit 68 ans!...
!
! Au-dela de cette duree, les entiers doivent etre 64 bits.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: IAMQ1
INTEGER(KIND=4) :: IAMQ2
INTEGER(KIND=4) :: IDELTA
INTEGER(KIND=4) :: IECJOURS
INTEGER(KIND=4) :: IH1
INTEGER(KIND=4) :: IH2
INTEGER(KIND=4) :: IM1
INTEGER(KIND=4) :: IM2
INTEGER(KIND=4) :: IRESTE
INTEGER(KIND=4) :: IS1
INTEGER(KIND=4) :: IS2
INTEGER(KIND=4) :: ISEC
INTEGER(KIND=4) :: KEC
character*(*) cd1,cd2
!
! -------------------------------------------------
! On lit les dates sur des entiers.
! -------------------------------------------------
!
read(cd1,fmt='(i8,3i2)') iamq1,ih1,im1,is1
!
! -------------------------------------------------
! Calculs d'ecarts et de leur partition
! en multiples de 86400 et sous-multiples.
! -------------------------------------------------
!
isec=ih1*3600+im1*60+is1 ! nombre de secondes ecoulees depuis cd10h.
idelta=kec+isec ! nombre de secondes entre cd10h et cd2.
ireste=modulo(idelta,86400) ! nombre de secondes entre cd20h et cd2.
iecjours=(idelta-ireste)/86400 ! nombre de jours entre cd10h et cd20h.
!
! -------------------------------------------------
! Date d'arrivee au jour pres.
! -------------------------------------------------
!
call daplus(iamq1,1,iecjours,iamq2)
!
! -------------------------------------------------
! Date finale a la seconde pres.
! -------------------------------------------------
!
ih2=ireste/3600
ireste=ireste-3600*ih2
im2=ireste/60
ireste=ireste-60*im2
is2=ireste
write(cd2,fmt='(i8,3i2.2)') iamq2,ih2,im2,is2
end
subroutine dj_vers_amqhms(pdj,kaaaa,kmm,kqq,khh,kmn,ps)
! --------------------------------------------------------------------------
! **** *dj_vers_amqhms*
! --------------------------------------------------------------------------
! Auteur:
! -------
! 1999-08-17, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! pdj date julienne associée à la date grégorienne UTC d'entrée
! En sortie:
! kaaaa année (4 chiffres!)
! kmm   mois
! kqq   quantième du mois
! khh   heure
! kmn   minute
! ps    seconde
! --------------------------------------------------------------------------
!
!-------------------------------------------------
! Nombre de jours entre le 1er janvier 2000 à 0 UTC
! et la date julienne courante.
!-------------------------------------------------
!


IMPLICIT NONE
INTEGER(KIND=4) :: IDATE1
INTEGER(KIND=4) :: IDATE2
INTEGER(KIND=4) :: IECART
INTEGER(KIND=4) :: KAAAA
INTEGER(KIND=4) :: KHH
INTEGER(KIND=4) :: KMM
INTEGER(KIND=4) :: KMN
INTEGER(KIND=4) :: KNOUV
INTEGER(KIND=4) :: KQQ
REAL(KIND=8) :: PDJ
REAL(KIND=8) :: PS
REAL(KIND=8) :: ZECART
REAL(KIND=8) :: ZFRAC
zecart=pdj-2451544.5
!
!-------------------------------------------------
! Nombre entier de jours.
!-------------------------------------------------
!
zfrac=modulo(zecart, 1._8 )
iecart=nint(zecart-zfrac)
!
!-------------------------------------------------
! Date grégorienne associée.
!-------------------------------------------------
!
idate1=20000101
call daplusj(idate1,iecart,idate2)
kqq=mod(idate2,100)
knouv=idate2/100
kmm=mod(knouv,100)
kaaaa=knouv/100
!
!-------------------------------------------------
! Calcul de des heure, minute et seconde.
!-------------------------------------------------
!
zfrac=(zecart-real(iecart))*24.
khh=int(zfrac)
zfrac=(zfrac-real(khh))*60.
kmn=int(zfrac)
ps=(zfrac-real(kmn))*60.
end
subroutine dj_vers_amqhmsree(pdj,pgrer)
! --------------------------------------------------------------------------
! **** **
! --------------------------------------------------------------------------
! Auteur:
! -------
! 2002-11, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! pdj date julienne
! En sortie:
! pgrer date en clair au format AAAAMMQQ.HHMMSS
! --------------------------------------------------------------------------
!

IMPLICIT NONE
REAL(KIND=8), intent(in) :: PDJ
REAL(KIND=8), intent(out) :: pgrer
REAL(KIND=8) :: ZS
INTEGER(KIND=4) :: iaaaa,imm,iqq,ihh,imn
!
!-------------------------------------------------
! Conversion grégorien julien; cible 5 entiers et un réel.
!-------------------------------------------------
!
call dj_vers_amqhms(pdj,iaaaa,imm,iqq,ihh,imn,zs)
!
!-------------------------------------------------
! On passe de ces 5 entiers et un réel à un seul réel.
!-------------------------------------------------
!
pgrer=real(iaaaa)*10000.+real(imm)*100. &
& + real(iqq)+real(ihh)/100. &
& + real(imn)/10000.+zs/1.E+06
end
subroutine amqhmsree_vers_dj(pgrer,pdj)
! --------------------------------------------------------------------------
! **** **
! --------------------------------------------------------------------------
! Auteur:
! -------
! 2002-11, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! pgrer date en clair au format AAAAMMQQ.HHMMSS
! En sortie:
! pdj date julienne associée à la date grégorienne
! --------------------------------------------------------------------------
!

IMPLICIT NONE
REAL(KIND=8), intent(out) :: PDJ
REAL(KIND=8), intent(in) :: pgrer
REAL(KIND=8) :: ZS,zloc
INTEGER(KIND=4) :: iaaaa,imm,iqq,ihh,imn,iloc
!
!-------------------------------------------------
! On passe de cette date grégorienne donnée
! comme un seul réel à 5 entiers et un réel.
!-------------------------------------------------
!
iloc=int(pgrer)
iqq=mod(iloc,100)
iloc=iloc/100
imm=mod(iloc,100)
iaaaa=iloc/100

iloc=nint((pgrer-real(int(pgrer)))*1.E+06)
zs=real(mod(iloc,100))
iloc=iloc/100
imn=mod(iloc,100)
ihh=iloc/100
!
!-------------------------------------------------
! Conversion grégorien julien; cible 5 entiers et un réel.
!-------------------------------------------------
!
call amqhms_vers_dj(iaaaa,imm,iqq,ihh,imn,zs,pdj)
end
subroutine ecartd(kdat1,kdat2,kopt,kgre)
! --------------------------------------------------------------------------
! **** *ECART* Ecart entre deux dates.
! --------------------------------------------------------------------------
! Auteur:
! -------
! 97-01-09, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree: kopt option de precision sur les dates:
! 1 : au jour pres
! 2 : a l'heure pres
! 3 : a la minute pres
! 4 : a la seconde pres
! - si kopt=1 : kdat au format AAAAMMQQ
! - si kopt=2 : kdat au format AAAAMMQQHH
! - si kopt=3 : kdat au format AAAAMMQQHHMM
! - si kopt=4 : kdat au format AAAAMMQQHHMMSS
! kdat1 et kdat2 dates au format ci-dessus.
! En sortie:
! - si kopt=1 : kgre nombre de jours    entre kdat1 et kdat2
! - si kopt=2 : kgre nombre d'heures    entre kdat1 et kdat2
! - si kopt=3 : kgre nombre de minutes  entre kdat1 et kdat2
! - si kopt=4 : kgre nombre de secondes entre kdat1 et kdat2
! --------------------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! 1. Vous compilez les entiers sur 32 bits:
! Vous devez alors vous limiter a kopt <= 2.
! L'ecart entre les deux dates doit etre inferieur a
! - 2**31 heures si kopt=2
! - 2**31 jours si kopt=1
! 2. Vous compilez les entiers sur 64 bits:
! Vous pouvez utiliser toutes les valeurs de kopt.
! L'ecart entre les deux dates doit etre inferieur a
! - 2**63 secondes si kopt=4
! - 2**63 minutes si kopt=3
! - 2**63 heures si kopt=2
! - 2**63 jours si kopt=1
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: IGRE1
INTEGER(KIND=4) :: IGRE2
INTEGER(KIND=4) :: KDAT1
INTEGER(KIND=4) :: KDAT2
INTEGER(KIND=4) :: KGRE
INTEGER(KIND=4) :: KOPT
call gregod(kdat1,kopt,igre1)
call gregod(kdat2,kopt,igre2)
kgre=igre2-igre1
end
subroutine ecartdj(k1,k2,kec)
! --------------------------------------------------------------
! **** *ecartdj* Ecart en jours entre deux dates.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   97-11, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! k1 date de depart au format AAAAMMQQ.
! k2 date d'arrivee au format AAAAMMQQ.
! En sortie:
! kec: nombre de jours entre les deux dates.
! --------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! Cette routine est utilisable avec des entiers 32 bits,
! si l'ecart entre les deux dates est inferieur a 2**31 jours,
! soit 5879489 ans!...
!
! Au-dela de cette duree, les entiers doivent etre 64 bits.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! -------------------------------------------------
! Ecart entre les deux dates au jour pres.
! -------------------------------------------------
!


IMPLICIT NONE
INTEGER(KIND=4) :: K1
INTEGER(KIND=4) :: K2
INTEGER(KIND=4) :: KEC
call ecartd(k1,k2,1,kec)
end
subroutine ecartds(cd1,cd2,kec)
! --------------------------------------------------------------
! **** *ecartds* Ecart en secondes entre deux dates.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   97-11, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! cd1 date de depart au format 'AAAAMMQQHHNNSS'.
! cd2 date d'arrivee au format 'AAAAMMQQHHNNSS'.
! En sortie:
! kec: nombre de secondes entre les deux dates.
! --------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! Cette routine est utilisable avec des entiers 32 bits,
! si l'ecart entre les deux dates est inferieur a 2**31 secondes,
! soit 68 ans!...
!
! Au-dela de cette duree, les entiers doivent etre 64 bits.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: IAMQ1
INTEGER(KIND=4) :: IAMQ2
INTEGER(KIND=4) :: IH1
INTEGER(KIND=4) :: IH2
INTEGER(KIND=4) :: IM1
INTEGER(KIND=4) :: IM2
INTEGER(KIND=4) :: IS1
INTEGER(KIND=4) :: IS2
INTEGER(KIND=4) :: KEC
INTEGER(KIND=4) :: KECQ
character*(*) cd1,cd2
!
! -------------------------------------------------
! On lit les dates sur des entiers.
! -------------------------------------------------
!
read(cd1,fmt='(i8,3i2)') iamq1,ih1,im1,is1
read(cd2,fmt='(i8,3i2)') iamq2,ih2,im2,is2
!
! -------------------------------------------------
! Ecart entre les deux dates au jour pres.
! -------------------------------------------------
!
call ecartd(iamq1,iamq2,1,kecq)
!
! -------------------------------------------------
! Ecart en secondes.
! -------------------------------------------------
!
kec=kecq*86400+(ih2-ih1)*3600+(im2-im1)*60+is2-is1
end
subroutine gregod(kdat,kopt,kgre)
! --------------------------------------------------------------------------
! **** *GREGOD *  - Conversion Date > Ecart par rapport a une date fixe.
! --------------------------------------------------------------------------
! Auteur:
! -------
! 92-05-27, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree: kopt option de precision sur les dates:
! 1 : au jour pres
! 2 : a l'heure pres
! 3 : a la minute pres
! 4 : a la seconde pres
! - si kopt=1 : kdat au format AAAAMMQQ
! - si kopt=2 : kdat au format AAAAMMQQHH
! - si kopt=3 : kdat au format AAAAMMQQHHMM
! - si kopt=4 : kdat au format AAAAMMQQHHMMSS
! En sortie:
! - si kopt=1 : kgre nombre de jours    entre 19000101       et kdat
! - si kopt=2 : kgre nombre d'heures    entre 1900010100     et kdat
! - si kopt=3 : kgre nombre de minutes  entre 190001010000   et kdat
! - si kopt=4 : kgre nombre de secondes entre 19000101000000 et kdat
! --------------------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! 1. Vous compilez les entiers sur 32 bits:
! Vous devez alors vous limiter a kopt <= 2.
! 2. Vous compilez les entiers sur 64 bits:
! Vous pouvez utiliser toutes les valeurs de kopt.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: idebm(12)
INTEGER(KIND=4) :: I0
INTEGER(KIND=4) :: IA100
INTEGER(KIND=4) :: IA4
INTEGER(KIND=4) :: IA400
INTEGER(KIND=4) :: IAAAA
INTEGER(KIND=4) :: IBISSEXT
INTEGER(KIND=4) :: ICONV
INTEGER(KIND=4) :: IFRJOUR
INTEGER(KIND=4) :: II
INTEGER(KIND=4) :: II1
INTEGER(KIND=4) :: IJOURP
INTEGER(KIND=4) :: IMM
INTEGER(KIND=4) :: IN
INTEGER(KIND=4) :: IN1
INTEGER(KIND=4) :: IN2
INTEGER(KIND=4) :: IQQ
INTEGER(KIND=4) :: KDAT
INTEGER(KIND=4) :: KGRE
INTEGER(KIND=4) :: KOPT
data idebm/0,31,59,90,120,151,181,212,243,273,304,334/
!
! --------------------------------------------------------------------------
! **      1. Calcul du nb de jours separant ki2 du 1er janv 1900
!
! *       1.1 Extraction des quantieme, mois et annee
if(kopt == 1) then
  ! Date de type AAAAMMQQ
  iconv=1
  ifrjour=0
  ii=kdat
elseif(kopt == 2) then
  ! Date de type AAAAMMQQHH
  iconv=24
  ifrjour=mod(kdat,100)
  ii=kdat/100
elseif(kopt == 3) then
  ! Date de type AAAAMMQQHHMM
  iconv=1440
  ifrjour=mod(kdat,100)
  ii=kdat/100
  ifrjour=ifrjour+mod(ii,100)*60
  ii=ii/100
elseif(kopt == 4) then
  ! Date de type AAAAMMQQHHMMSS
  iconv=86400
  ifrjour=mod(kdat,100)
  ii=kdat/100
  ifrjour=ifrjour+mod(ii,100)*60
  ii=ii/100
  ifrjour=ifrjour+mod(ii,100)*3600
  ii=ii/100
else
  ! Cas d'entree erronee de l'utilisateur.
  print*,'GREGODR/ERREUR: argument kopt errone!...'
  print*,kopt
  stop 'call abort'
endif
iqq=ii-(ii/100)*100
in=(ii-iqq)/100
imm=in-(in/100)*100
iaaaa=(in-imm)/100
! *       1.2 L'annee est-elle bissextile?
! Une annee est bissextile ssi elle est
! (mult de 4 et non mult de 100) ou (mult de 400)
iaaaa=iaaaa
ia400=400*(iaaaa/400)
ia100=100*(iaaaa/100)
ia4=4*(iaaaa/4)
if((iaaaa == ia400).or.((iaaaa == ia4).and.(iaaaa /= ia100)))then
  ibissext=1
else
  ibissext=0
endif
if ((ibissext == 1).and.(imm > 2)) then
  ijourp=1
else
  ijourp=0
endif
! *       1.3 Nombre de jours ecoules depuis le 1er janv
if(imm > 12) then
  print*,'GREGODR/ERREUR: mois errone.'
  print*,imm
  stop 'call abort'
endif
in2=idebm(imm)+ijourp+iqq-1
! *       1.4 Calcul du nb de jours separant les 1er janvier de ii et 1900
i0=1900
in2=in2+365*(iaaaa-i0)+int((iaaaa-1)/4)-int((i0-1)/4)&
&-int((iaaaa-1)/100)+int((i0-1)/100)&
&+int((iaaaa-1)/400)-int((i0-1)/400)
! --------------------------------------------------------------------------
! **      2. Calcul du nb de jours separant ii1 du 1er janv 1900
!
! *       2.1 Extraction des quantieme, mois et annee
ii1=19000101
ii=ii1
iqq=ii-(ii/100)*100
in=(ii-iqq)/100
imm=in-(in/100)*100
iaaaa=(in-imm)/100
! *       2.2 L'annee est-elle bissextile?
! Une annee est bissextile ssi elle est
! (mult de 4 et non mult de 100) ou (mult de 400)
iaaaa=iaaaa
ia400=400*(iaaaa/400)
ia100=100*(iaaaa/100)
ia4=4*(iaaaa/4)
if((iaaaa == ia400).or.((iaaaa == ia4).and.(iaaaa /= ia100)))then
  ibissext=1
else
  ibissext=0
endif
if ((ibissext == 1).and.(imm > 2)) then
  ijourp=1
else
  ijourp=0
endif
! *       2.3 Nombre de jours ecoules depuis le 1er janv
in1=idebm(imm)+ijourp+iqq-1
! *       2.4 Calcul du nb de jours separant les 1er janvier de ii et 1900
i0=1900
in1=in1+365*(iaaaa-i0)+int((iaaaa-1)/4)-int((i0-1)/4)&
&-int((iaaaa-1)/100)+int((i0-1)/100)&
&+int((iaaaa-1)/400)-int((i0-1)/400)
! --------------------------------------------------------------------------
! **      3. Difference in2-in1
kgre=(in2-in1)*iconv+ifrjour
end
subroutine gregoi(kgre,kopt,kdat)
! --------------------------------------------------------------------------
! **** *GREGOI *  - Conversion Ecart par rapport a une date fixe > Date.
! --------------------------------------------------------------------------
! Auteur:
! -------
! 92-05-27, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree: kopt option de precision sur les dates:
! 1 : au jour pres
! 2 : a l'heure pres
! 3 : a la minute pres
! 4 : a la seconde pres
! - si kopt=1 : kgre nombre de jours    entre 19000101       et kdat
! - si kopt=2 : kgre nombre d'heures    entre 1900010100     et kdat
! - si kopt=3 : kgre nombre de minutes  entre 190001010000   et kdat
! - si kopt=4 : kgre nombre de secondes entre 19000101000000 et kdat
! En sortie:
! - si kopt=1 : kdat au format AAAAMMQQ
! - si kopt=2 : kdat au format AAAAMMQQHH
! - si kopt=3 : kdat au format AAAAMMQQHHMM
! - si kopt=4 : kdat au format AAAAMMQQHHMMSS
! --------------------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! 1. Vous compilez les entiers sur 32 bits:
! Vous devez alors vous limiter a kopt <= 2.
! 2. Vous compilez les entiers sur 64 bits:
! Vous pouvez utiliser toutes les valeurs de kopt.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: ijours(12)
INTEGER(KIND=4) :: IA100
INTEGER(KIND=4) :: IA4
INTEGER(KIND=4) :: IA400
INTEGER(KIND=4) :: IAAAA
INTEGER(KIND=4) :: IBISSEXT
INTEGER(KIND=4) :: ICONV
INTEGER(KIND=4) :: IDAT
INTEGER(KIND=4) :: IEC
INTEGER(KIND=4) :: IECI
INTEGER(KIND=4) :: IGII2P
INTEGER(KIND=4) :: II2P
INTEGER(KIND=4) :: IMM
INTEGER(KIND=4) :: IMOD
INTEGER(KIND=4) :: IQQ
INTEGER(KIND=4) :: KDAT
INTEGER(KIND=4) :: KGRE
INTEGER(KIND=4) :: KOPT
REAL(KIND=8) :: ZARRDEC
data ijours/31,28,31,30,31,30,31,31,30,31,30,31/
! --------------------------------------------------------------------------
! **   On determine la date approximative d'arrivee en annees decimales
!
if(kopt == 1) then
  ! Date de type AAAAMMQQ
  iconv=1
elseif(kopt == 2) then
  ! Date de type AAAAMMQQHH
  iconv=24
elseif(kopt == 3) then
  ! Date de type AAAAMMQQHHMM
  iconv=1440
elseif(kopt == 4) then
  ! Date de type AAAAMMQQHHMMSS
  iconv=86400
else
  ! Cas d'entree erronee de l'utilisateur.
  print*,'GREGOI/ERREUR: argument kopt errone!...'
  print*,kopt
  stop 'call abort'
endif
zarrdec=1900.+(real(kgre)/real(iconv)-5.)/365.2425
! --------------------------------------------------------------------------
! **   On determine la date en clair ii2p associee a la date decimale
!
iaaaa=int(zarrdec)
zarrdec=12.*(zarrdec-real(iaaaa))
imm=int(zarrdec)+1
zarrdec=28.*(zarrdec-real(imm-1))
iqq=int(zarrdec)+1
ii2p=iqq+imm*100+iaaaa*10000
! --------------------------------------------------------------------------
! **   On calcule le nombre de jours separant 19000101 de ii2p
!
call gregod(ii2p,1,igii2p)
imod=mod(kgre,iconv)
if(imod < 0) imod=imod+iconv
iec=(kgre-imod)/iconv-igii2p
! --------------------------------------------------------------------------
! **   On avance de iec jours par rapport a ii2p
!
! *       L'annee est-elle bissextile?
! Une annee est bissextile ssi elle est
! (mult de 4 et non mult de 100) ou (mult de 400)
iaaaa=iaaaa
ia400=400*(iaaaa/400)
ia100=100*(iaaaa/100)
ia4=4*(iaaaa/4)
if((iaaaa == ia400).or.((iaaaa == ia4).and.(iaaaa /= ia100)))then
  ibissext=1
else
  ibissext=0
endif
! Si oui, 29 jours en fevrier
if(ibissext == 1) ijours(2)=29
! *       Boucle sur les jours
do ieci=1,iec
  iqq=iqq+1
  if(iqq > ijours(imm)) then
    iqq=1
    imm=imm+1
  endif
  if(imm > 12) then
    imm=1
    iaaaa=iaaaa+1
  endif
enddo
! --------------------------------------------------------------------------
! **   On met en forme la date finale
!
idat=iqq+imm*100+iaaaa*10000
if(kopt == 2) then
  imod=mod(kgre,iconv)
  if(imod < 0) imod=imod+iconv
  idat=idat*100+imod
elseif(kopt == 3) then
  imod=mod(kgre,iconv)
  if(imod < 0) imod=imod+iconv
  idat=idat*100+imod/60
  imod=mod(imod,60)
  idat=idat*100+imod
elseif(kopt == 4) then
  imod=mod(kgre,iconv)
  if(imod < 0) imod=imod+iconv
  idat=idat*100+imod/3600
  imod=mod(imod,3600)
  idat=idat*100+imod/60
  imod=mod(imod,60)
  idat=idat*100+imod
endif
kdat=idat
end
function ijoursem(kdat)
! --------------------------------------------------------------------------
! **** *IJOURSEM* Jour de la semaine de la date d'entree.
! --------------------------------------------------------------------------
! Auteur:
! -------
! 94-10-31, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! kdat1 au format AAAAMMQQ
! En sortie:
! ijour=0 si dimanche, 1 lundi, ..., 6 samedi.
! --------------------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: IDATDIM
INTEGER(KIND=4) :: IECART
INTEGER(KIND=4) :: IGRE
INTEGER(KIND=4) :: IGREDIM
INTEGER(KIND=4) :: KDAT
INTEGER(KIND=4) :: ijoursem
call gregod(kdat,1,igre)
idatdim=19941030 ! cette date etait un dimanche.
call gregod(idatdim,1,igredim)
iecart=igre-igredim
ijoursem=modulo(iecart,7)
end
subroutine qqmmaa(kdatd,cdresd)
! --------------------------------------------------------------------------
! **** *QQMMAA *  - Conversion d'un entier type AAAAQQMM vers une date en clair.
! --------------------------------------------------------------------------
! Auteur:
! -------
! 92-05-27, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: IAN
INTEGER(KIND=4) :: IGRE
INTEGER(KIND=4) :: ILOC
INTEGER(KIND=4) :: IMM
INTEGER(KIND=4) :: IQQ
INTEGER(KIND=4) :: KDATD
character*(*) cdresd
character*03 cljour
iqq=mod(kdatd,100)
iloc=kdatd/100
imm=mod(iloc,100)
ian=iloc/100
call gregod(kdatd,1,igre)
igre=mod(igre,7)
if(igre == 0) then
  cljour='Lun'
elseif(igre == 1) then
  cljour='Mar'
elseif(igre == 2) then
  cljour='Mer'
elseif(igre == 3) then
  cljour='Jeu'
elseif(igre == 4) then
  cljour='Ven'
elseif(igre == 5) then
  cljour='Sam'
elseif(igre == 6) then
  cljour='Dim'
endif
write(cdresd,fmt='(a3,a1,i2,a1,i2.2,a1,i4.4)')&
&cljour,' ',iqq,'.',imm,'.',ian
end
subroutine quant(kdate,kquant)
! --------------------------------------------------------------
! **** *quant* Quantieme de l'annee d'une date donnee.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   1999-04, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! kdate date au format AAAAMMQQ.
! En sortie:
! quantieme de l'annee (1 le 1er janvier, 32 le 1er fevrier, etc...)
! --------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: IBASE
INTEGER(KIND=4) :: IEC
INTEGER(KIND=4) :: KDATE
INTEGER(KIND=4) :: KQUANT
ibase=10000*(kdate/10000)+0101 ! 1er janvier de l'annee courante.
call ecartdj(ibase,kdate,iec)
kquant=iec+1
end
