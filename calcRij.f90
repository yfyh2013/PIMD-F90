      subroutine calcRij(Rij, R2, box)
         implicit none
         double precision, dimension(3) :: Rij, box
         double precision :: R2
         double precision, dimension(3) :: R1
         integer :: ii

            do ii=1, 3
               Rij(ii) = Rij(ii) + box (ii)
               do while (Rij(ii) .gt. 0.5d0*box (ii))
                  Rij(ii) = Rij(ii) - box (ii)
               end do
               do while (Rij(ii) .lt. -0.5d0*box (ii))
                  Rij(ii) = Rij(ii) + box (ii)
               end do
            enddo
          R2 = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
      end subroutine calcRij

