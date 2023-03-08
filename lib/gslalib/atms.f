      SUBROUTINE slATMS (RT, TT, DNT, GAMAL, R, DN, RDNDR)
*+
*     - - - - -
*      A T M S
*     - - - - -
*
*  Internal routine used by REFRO
*
*  Refractive index and derivative with respect to height for the
*  stratosphere.
*
*  Given:
*    RT      d    height of tropopause from centre of the Earth (metre)
*    TT      d    temperature at the tropopause (K)
*    DNT     d    refractive index at the tropopause
*    GAMAL   d    constant of the atmospheric model = G*MD/R
*    R       d    current distance from the centre of the Earth (metre)
*
*  Returned:
*    DN      d    refractive index at R
*    RDNDR   d    R * rate the refractive index is changing at R
*
*  Last revision:   26 December 2004
*
*-----------------------------------------------------------------------
*
*  ---------------------------------
*  Tpoint Software License Agreement
*  ---------------------------------
*
*  All software, both binary and source published by Tpoint Software
*  (hereafter, Software) is copyrighted by the the Company
*  (hereafter TPS), and ownership of all right, title and interest in
*  and to the Software remains with the TPS.  By using or copying the
*  Software, you (hereafter User) agrees to abide by the terms of this
*  Agreement.
*
*  Noncommercial Use
*
*  TPS grants to User a royalty-free, nonexclusive right to execute,
*  copy, modify and distribute both the binary and source code solely
*  for academic, research and other similar noncommercial uses, subject
*  to the following conditions:
*
*  1) User acknowledges that the Software is still under development and
*     that it is being supplied "as is," without any support services
*     from TPS.  Neither TPS nor the author makes any representations
*     or warranties, express or implied, including, without limitation,
*     any representations or warranties of the merchantability or
*     fitness for any particular purpose, or that the application of the
*     software, will not infringe on any patents or other proprietary
*     rights of others.
*
*  2) TPS shall not be held liable for direct, indirect, incidental or
*     consequential damages arising from any claim by User or any third
*     party with respect to uses allowed under this Agreement, or from
*     any use of the Software.
*
*  3) User agrees to indemnify fully and hold harmless TPS and/or the
*     author of the original work from and against any and all claims,
*     demands, suits, losses, damages, costs and expenses arising out of
*     User's use of the Software, including, without limitation, arising
*     out of User's modification of the Software.
*
*  4) User may modify the Software and distribute that modified work to
*     third parties provided that: (a) if posted separately, it clearly
*     acknowledges that it contains material copyrighted by TPS (b) no
*     charge is associated with such copies, and (d) User clearly
*     notifies secondary users that such modified work is not the
*     original Software.
*
*  5) User agrees that TPS, the authors of the original work and others
*     may enjoy a royalty-free, non-exclusive license to use, copy,
*     modify and redistribute these modifications to the Software made
*     by the User and distributed to third parties as a derivative work
*     under this agreement.
*
*  6) This agreement will terminate immediately upon User's breach of,
*     or non-compliance with, any of its terms.  User may be held liable
*     for any copyright infringement or the infringement of any other
*     proprietary rights in the Software that is caused or facilitated
*     by the User's failure to abide by the terms of this agreement.
*
*  7) This agreement will be construed and enforced in accordance with
*     English law applicable to contracts performed entirely within the
*     United Kingdom.  The parties irrevocably consent to the exclusive
*     jurisdiction of the English courts for all disputes concerning
*     this agreement.
*
*  Commercial Use
*
*  Any User wishing to make a commercial use of the Software must
*  contact TPS to arrange an appropriate license.  Commercial use
*  includes (1) integrating or incorporating all or part of the source
*  code and/or intellectual content not published elsewhere into a
*  product for sale or license by, or on behalf of, User to third
*  parties, or (2) distribution of the binary code or source code or
*  intellectual content not published elsewhere to third parties for
*  use with a commercial product sold or licensed by, or on behalf of,
*  User.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION RT, TT, DNT, GAMAL, R, DN, RDNDR

      DOUBLE PRECISION B, W


      B = GAMAL/TT
      W = (DNT-1D0)*EXP(-B*(R-RT))
      DN = 1D0+W
      RDNDR = -R*B*W

      END
