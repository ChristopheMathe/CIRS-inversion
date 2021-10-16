subroutine database_planet(planet_name,lat,g_1bar,r_1bar)

implicit none
double precision, parameter :: dpi = 1.745329252d-02
double precision :: sin2 = 0.0d0
double precision :: ar = 0.0d0
double precision :: planet_radii
double precision :: planet_gravity
double precision :: planet_flatness !flatness
double precision :: planet_grav_field_rotation !gravity field q(rotation)
double precision :: planet_grav_field_momentum !gravity field j2(momentum)
double precision :: lat
character(len=256) :: planet_name

double precision, intent(out) :: g_1bar
double precision, intent(out) :: r_1bar

select case(planet_name)

    case('TIT')
        planet_radii = 2575.
        planet_gravity = 135.!.4
        planet_flatness = 0.
        planet_grav_field_rotation = 0. !q
        planet_grav_field_momentum = 0. !j2
    case default
        print *, planet_name, "is not referred in database_planet.f90"
        stop
end select

!**             Determination of the planet and gravity at surface
sin2=dsin(dpi*lat)**2
ar=dsqrt(1.+sin2*planet_flatness*(2.-planet_flatness)/(1.-planet_flatness)**2)
g_1bar = planet_gravity*ar*ar*(1.-1.5*planet_grav_field_momentum*(3.*sin2-1.)*ar*ar&
         -planet_grav_field_rotation*(1.-sin2)/(ar*ar*ar))/(1.+1.5*planet_grav_field_momentum-planet_grav_field_rotation)
r_1bar=planet_radii/ar

end subroutine database_planet
