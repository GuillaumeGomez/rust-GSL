//
// A rust binding for the GSL library by Guillaume Gomez (guillaume1.gomez@gmail.com)
//

// The following program demonstrates the use of the physical constants in a calculation. In this case, the goal is to calculate the range
// of light-travel times from Earth to Mars.
//
// The required data is the average distance of each planet from the Sun in astronomical units (the eccentricities and inclinations of the
// orbits will be neglected for the purposes of this calculation). The average radius of the orbit of Mars is 1.52 astronomical units, and
// for the orbit of Earth it is 1 astronomical unit (by definition). These values are combined with the MKSA values of the constants for
// the speed of light and the length of an astronomical unit to produce a result for the shortest and longest light-travel times in
// seconds. The figures are converted into minutes before being displayed.

extern crate rgsl;

fn main() {
    let c = rgsl::physical_constant::MKSA_SPEED_OF_LIGHT;
    let au = rgsl::physical_constant::MKSA_ASTRONOMICAL_UNIT;
    let minutes = rgsl::physical_constant::MKSA_MINUTE;

    /* distance stored in meters */
    let r_earth = 1f64 * au;
    let r_mars = 1.52f64 * au;

    let t_min = (r_mars - r_earth) / c;
    let t_max = (r_mars + r_earth) / c;

    println!("light travel time from Earth to Mars:");
    println!("minimum = {:.1} minutes", t_min / minutes);
    println!("maximum = {:.1} minutes", t_max / minutes);
}
