package com.hwaipy.predict

import java.io.File
import scala.collection.mutable.ArrayBuffer

class Predict {
  private val deg2rad = 1.745329251994330E-2
  /* Degrees to radians */
  private val pio2 = Math.PI / 2
  /* Pi/2 */
  private val x3pio2 = pio2 * 3
  /* 3*Pi/2 */
  private var twopi = Math.PI * 2
  /* 2*Pi  */
  private val e6a = 1.0E-6
  private val tothrd = 6.6666666666666666E-1
  /* 2/3 */
  private val xj2 = 1.0826158E-3
  /* J2 Harmonic (WGS '72) */
  private val xj3 = -2.53881E-6
  /* J3 Harmonic (WGS '72) */
  private val xj4 = -1.65597E-6
  /* J4 Harmonic (WGS '72) */
  private val xke = 7.43669161E-2
  private val xkmper = 6.378137E3
  /* WGS 84 Earth radius km */
  private val xmnpda = 1.44E3
  /* Minutes per day */
  private val ae = 1.0
  private val ck2 = 5.413079E-4
  private val ck4 = 6.209887E-7
  private val f = 3.35281066474748E-3
  /* Flattening factor */
  private val ge = 3.986008E5
  /* Earth gravitational constant (WGS '72) */
  private val s = 1.012229
  private val qoms2t = 1.880279E-09
  private val secday = 8.6400E4
  /* Seconds per day */
  private val omega_E = 1.00273790934
  /* Earth rotations/siderial day */
  private val omega_ER = 6.3003879
  /* Earth rotations, rads/siderial day */
  private val zns = 1.19459E-5
  private val c1ss = 2.9864797E-6
  private val zes = 1.675E-2
  private val znl = 1.5835218E-4
  private val c1l = 4.7968065E-7
  private val zel = 5.490E-2
  private val zcosis = 9.1744867E-1
  private val zsinis = 3.9785416E-1
  private val zsings = -9.8088458E-1
  private val zcosgs = 1.945905E-1
  private val zcoshs = 1
  private val zsinhs = 0
  private val q22 = 1.7891679E-6
  private val q31 = 2.1460748E-6
  private val q33 = 2.2123015E-7
  private val g22 = 5.7686396
  private val g32 = 9.5240898E-1
  private val g44 = 1.8014998
  private val g52 = 1.0508330
  private val g54 = 4.4108898
  private val root22 = 1.7891679E-6
  private val root32 = 3.7393792E-7
  private val root44 = 7.3636953E-9
  private val root52 = 1.1428639E-7
  private val root54 = 2.1765803E-9
  private val thdt = 4.3752691E-3
  private val rho = 1.5696615E-1
  private val mfactor = 7.292115E-5
  private val sr = 6.96000E5
  /* Solar radius - km (IAU 76) */
  private val AU = 1.49597870691E8 /* Astronomical unit - km (IAU 76) */

  /* Entry points of Deep() */

  private var dpinit = 1
  /* Deep-space initialization code */
  private var dpsec = 2
  /* Deep-space secular code        */
  private var dpper = 3 /* Deep-space periodic code       */

  /* Flow control flag definitions */

  private var ALL_FLAGS = -1
  private var SGP_INITIALIZED_FLAG = 0x000001
  /* not used */
  private var SGP4_INITIALIZED_FLAG = 0x000002
  private var SDP4_INITIALIZED_FLAG = 0x000004
  private var SGP8_INITIALIZED_FLAG = 0x000008
  private var SDP8_INITIALIZED_FLAG = 0x000010
  private var SIMPLE_FLAG = 0x000020
  private var DEEP_SPACE_EPHEM_FLAG = 0x000040
  private var LUNAR_TERMS_DONE_FLAG = 0x000080
  private var NEW_EPHEMERIS_FLAG = 0x000100
  private var DO_LOOP_FLAG = 0x000200
  private var RESONANCE_FLAG = 0x000400
  private var SYNCHRONOUS_FLAG = 0x000800
  private var EPOCH_RESTART_FLAG = 0x001000
  private var VISIBLE_FLAG = 0x002000
  private var SAT_ECLIPSED_FLAG = 0x004000

  class Sat {
    var line1 = ""
    var line2 = ""
    var name = ""
    var catnum = 0L
    var setnum = 0L
    var designator = ""
    var year = 0
    var refepoch = .0
    var incl = .0
    var raan = .0
    var eccn = .0
    var argper = .0
    var meanan = .0
    var meanmo = .0
    var drag = .0
    var nddot6 = .0
    var bstar = .0
    var orbitnum = 0L
  }

  var sat = new Sat

  class Qth {
    var callsign = new Array[Char](17)
    var stnlat = .0
    var stnlong = .0
    var stnalt = 0
  }

  var qth = new Qth

  /* Global variables for sharing data among functions... */
  var tsince = .0
  var jul_epoch = .0
  var jul_utc = .0
  var eclipse_depth = 0
  var sat_azi = .0
  var sat_ele = .0
  var sat_range = .0
  var sat_range_rate = .0
  var sat_lat = .0
  var sat_lon = .0
  var sat_alt = .0
  var sat_vel = .0
  var phase = .0
  var sun_azi = .0
  var sun_ele = .0
  var daynum = .0
  var fm = .0
  var fk = .0
  var age = .0
  var aostime = .0
  var lostime = .0
  var ax = .0
  var ay = .0
  var az = .0
  var rx = .0
  var ry = .0
  var rz = .0
  var squint = .0
  var alat = .0
  var alon = .0
  var sun_ra = .0
  var sun_dec = .0
  var sun_lat = .0
  var sun_lon = .0
  var sun_range = .0
  var sun_range_rate = .0
  var moon_az = .0
  var moon_el = .0
  var moon_dx = .0
  var moon_ra = .0
  var moon_dec = .0
  var moon_gha = .0
  var temp = new Array[Char](80)
  var output = new Array[Char](25)
  var serial_port = new Array[Char](15)
  var netport = new Array[Char](6)
  var ephem = ""
  var reload_tle = 0
  var once_per_second = 0
  var sat_sun_status = 0
  var findsun = 0
  var calc_squint = 0
  var database = 0
  var io_lat = 'N'
  var io_lon = 'W'
  var antfd = 0
  var iaz = 0
  var iel = 0
  var ma256 = 0
  var isplat = 0
  var isplong = 0
  var socket_flag = 0
  var Flags = 0
  var fel = 0
  var faz = 0
  var fdistance = 0
  var rv = 0L
  var irk = 0L
  var valval = new Array[Int](256)

  /* The following variables are used by the socket server.  They
 are updated in the MultiTrack() and SingleTrack() functions. */ var visibility_array = new Array[Char](24)
  var tracking_mode = new Array[Char](30)

  var az_array = new Array[Double](24)
  var el_array = new Array[Double](24)
  var long_array = new Array[Double](24)
  var lat_array = new Array[Double](24)
  var footprint_array = new Array[Double](24)
  var altitude_array = new Array[Double](24)
  var velocity_array = new Array[Double](24)
  var eclipse_depth_array = new Array[Double](24)
  var phase_array = new Array[Double](24)
  var squint_array = new Array[Double](24)

  var doppler = new Array[Double](24)
  var nextevent = new Array[Double](24)

  var aos_array = new Array[Long](24)
  var orbitnum_array = new Array[Long](24)

  var portbase = 0


  /**
    * Type definitions *
    */
  /* Two-line-element satellite orbital data
     structure used directly by the SGP4/SDP4 code. */ class tle_t {
    var epoch = .0
    var xndt2o = .0
    var xndd6o = .0
    var bstar = .0
    var xincl = .0
    var xnodeo = .0
    var eo = .0
    var omegao = .0
    var xmo = .0
    var xno = .0
    var catnr = 0L
    var elset = 0L
    var revnum = 0L
    var sat_name = ""
    var idesg = ""
  }

  /* Geodetic position structure used by SGP4/SDP4 code. */ class geodetic_t {
    var lat = .0
    var lon = .0
    var alt = .0
    var theta = .0
  }

  /* General three-dimensional vector structure used by SGP4/SDP4 code. */ class vector_t {
    var x = .0
    var y = .0
    var z = .0
    var w = .0
  }


  /* Common arguments between deep-space functions used by SGP4/SDP4 code. */
  class deep_arg_t {
    /* Used by dpinit part of Deep() */ var eosq = .0
    var sinio = .0
    var cosio = .0
    var betao = .0
    var aodp = .0
    var theta2 = .0
    var sing = .0
    var cosg = .0
    var betao2 = .0
    var xmdot = .0
    var omgdot = .0
    var xnodot = .0
    var xnodp = .0
    /* Used by dpsec and dpper parts of Deep() */
    var xll = .0
    var omgadf = .0
    var xnode = .0
    var em = .0
    var xinc = .0
    var xn = .0
    var t = .0
    /* Used by thetg and Deep() */
    var ds50 = .0
  }

  /* Global structure used by SGP4/SDP4 code. */
  var obs_geodetic = new geodetic_t

  /* Two-line Orbital Elements for the satellite used by SGP4/SDP4 code. */
  var tle = new tle_t

  /* Functions for testing and setting/clearing flags used in SGP4/SDP4 code */
  def isFlagSet(flag: Int): Boolean = return (Flags & flag) != 0

  def isFlagClear(flag: Int): Boolean = {
    return (~(Flags) & flag) != 0
  }

  def SetFlag(flag: Int): Unit = {
    Flags |= flag
  }

  def ClearFlag(flag: Int): Unit = {
    Flags &= ~(flag)
  }

  /* Remaining SGP4/SDP4 code follows... */
  def Sign(arg: Double): Int = {
    /* Returns sign of a double */ if (arg > 0) {
      return 1
    }
    else {
      if (arg < 0) {
        return -(1)
      }
      else {
        return 0
      }
    }
  }

  def Sqr(arg: Double): Double = /* Returns square of a double */ arg * arg

  def Cube(arg: Double): Double = /* Returns cube of a double */ arg * arg * arg

  def Radians(arg: Double): Double = /* Returns angle in radians from argument in degrees */ arg * deg2rad

  def Degrees(arg: Double): Double = /* Returns angle in degrees from argument in radians */ arg / deg2rad

  def ArcSin(arg: Double): Double = /* Returns the arcsine of the argument */ if (Math.abs(arg) >= 1.0) Sign(arg) * pio2
  else Math.atan(arg / Math.sqrt(1.0 - arg * arg))

  def ArcCos(arg: Double): Double = /* Returns arccosine of argument */ pio2 - ArcSin(arg)

  def Magnitude(v: Predict#vector_t): Unit = {
    /* Calculates scalar magnitude of a vector_t argument */ v.w = Math.sqrt(Sqr(v.x) + Sqr(v.y) + Sqr(v.z))
  }

  def Vec_Add(v1: Predict#vector_t, v2: Predict#vector_t, v3: Predict#vector_t): Unit = {
    /* Adds vectors v1 and v2 together to produce v3 */ v3.x = v1.x + v2.x
    v3.y = v1.y + v2.y
    v3.z = v1.z + v2.z
    Magnitude(v3)
  }

  def Vec_Sub(v1: Predict#vector_t, v2: Predict#vector_t, v3: Predict#vector_t): Unit = {
    /* Subtracts vector v2 from v1 to produce v3 */ v3.x = v1.x - v2.x
    v3.y = v1.y - v2.y
    v3.z = v1.z - v2.z
    Magnitude(v3)
  }

  def Scalar_Multiply(k: Double, v1: Predict#vector_t, v2: Predict#vector_t): Unit = {
    /* Multiplies the vector v1 by the scalar k to produce the vector v2 */ v2.x = k * v1.x
    v2.y = k * v1.y
    v2.z = k * v1.z
    v2.w = Math.abs(k) * v1.w
  }

  def Scale_Vector(k: Double, v: Predict#vector_t): Unit = {
    /* Multiplies the vector v1 by the scalar k */ v.x *= k
    v.y *= k
    v.z *= k
    Magnitude(v)
  }

  def Dot(v1: Predict#vector_t, v2: Predict#vector_t): Double = /* Returns the dot product of two vectors */ v1.x * v2.x + v1.y * v2.y + v1.z * v2.z

  def Angle(v1: Predict#vector_t, v2: Predict#vector_t): Double = {
    /* Calculates the angle between vectors v1 and v2 */ Magnitude(v1)
    Magnitude(v2)
    ArcCos(Dot(v1, v2) / (v1.w * v2.w))
  }

  def Cross(v1: Predict#vector_t, v2: Predict#vector_t, v3: Predict#vector_t): Unit = {
    /* Produces cross product of v1 and v2, and returns in v3 */ v3.x = v1.y * v2.z - v1.z * v2.y
    v3.y = v1.z * v2.x - v1.x * v2.z
    v3.z = v1.x * v2.y - v1.y * v2.x
    Magnitude(v3)
  }

  def Normalize(v: Predict#vector_t): Unit = {
    /* Normalizes a vector */ v.x /= v.w
    v.y /= v.w
    v.z /= v.w
  }

  def AcTan(sinx: Double, cosx: Double): Double = /* Four-quadrant arctan function */ if (cosx == 0.0) if (sinx > 0.0) pio2
  else x3pio2
  else if (cosx > 0.0) if (sinx > 0.0) Math.atan(sinx / cosx)
  else twopi + Math.atan(sinx / cosx)
  else Math.PI + Math.atan(sinx / cosx)

  def FMod2p(x: Double): Double = {
    /* Returns mod 2PI of argument */ var i = 0
    var ret_val = .0
    ret_val = x
    i = (ret_val / twopi).asInstanceOf[Int]
    ret_val -= i * twopi
    if (ret_val < 0.0) ret_val += twopi
    ret_val
  }

  def Modulus(arg1: Double, arg2: Double): Double = {
    /* Returns arg1 mod arg2 */ var i = 0
    var ret_val = .0
    ret_val = arg1
    i = (ret_val / arg2).toInt
    ret_val -= i * arg2
    if (ret_val < 0.0) ret_val += arg2
    ret_val
  }

  def Frac(arg: Double): Double = /* Returns fractional part of double argument */ arg - Math.floor(arg)

  def Round(arg: Double): Int = /* Returns argument rounded up to nearest integer */ Math.floor(arg + 0.5).toInt

  def Int(arg: Double): Double = /* Returns the floor integer of a double arguement, as double */ Math.floor(arg)

  def Convert_Sat_State(pos: Predict#vector_t, vel: Predict#vector_t): Unit = {
    /* Converts the satellite's position and velocity  *//* vectors from normalized values to km and km/sec */ Scale_Vector(xkmper, pos)
    Scale_Vector(xkmper * xmnpda / secday, vel)
  }

  def Julian_Date_of_Year(yearPara: Double): Double = {
    /* The function Julian_Date_of_Year calculates the Julian Date  *//* of Day 0.0 of {year}. This function is used to calculate the *//* Julian Date of any date by using Julian_Date_of_Year, DOY,   *//* and Fraction_of_Day. *//* Astronomical Formulae for Calculators, Jean Meeus, *//* pages 23-25. Calculate Julian Date of 0.0 Jan year */ var A = 0L
    var B = 0L
    var i = 0L
    var jdoy = .0
    var year = yearPara - 1
    i = (year / 100).toLong
    A = i
    i = A / 4
    B = 2 - A + i
    i = (365.25 * year).toLong
    i += (30.6001 * 14).toLong
    jdoy = i + 1720994.5 + B
    jdoy
  }

  def modf_Integer(x: Double): Double = x.toInt

  def modf_Frac(x: Double): Double = x - x.toInt

  def Julian_Date_of_Epoch(epoch: Double): Double = {
    /* The function Julian_Date_of_Epoch returns the Julian Date of     *//* an epoch specified in the format used in the NORAD two-line      *//* element sets. It has been modified to support dates beyond       *//* the year 1999 assuming that two-digit years in the range 00-56   *//* correspond to 2000-2056. Until the two-line element set format   *//* is changed, it is only valid for dates through 2056 December 31. */ var year = .0
    var day = .0
    /* Modification to support Y2K *//* Valid 1957 through 2056     */ year = modf_Integer(epoch * 1E-3)
    day = modf_Frac(epoch * 1E-3) * 1E3
    if (year < 57) year = year + 2000
    else year = year + 1900
    Julian_Date_of_Year(year) + day
  }

  def DOY(yr: Int, mo: Int, dy: Int): Int = {
    /* The function DOY calculates the day of the year for the specified *//* date. The calculation uses the rules for the Gregorian calendar   *//* and is valid from the inception of that calendar system.          */ val days = Array(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    var i = 0
    var day = 0
    day = 0
    i = 0
    while ( {
      i < mo - 1
    }) {
      day += days(i)

      {
        i += 1;
        i - 1
      }
    }
    day = day + dy
    /* Leap year correction */ if ((yr % 4 == 0) && ((yr % 100 != 0) || (yr % 400 == 0)) && (mo > 2)) day += 1
    day
  }

  def Fraction_of_Day(hr: Int, mi: Int, se: Double): Double = {
    /* Fraction_of_Day calculates the fraction of *//* a day passed at the specified input time.  */ var dhr = .0
    var dmi = .0
    dhr = hr.toDouble
    dmi = mi.toDouble
    (dhr + (dmi + se / 60.0) / 60.0) / 24.0
  }

  class tm {
    val tm_sec = 0
    /* 秒 – 取值区间为[0,59] */
    val tm_min = 0
    /* 分 - 取值区间为[0,59] */
    val tm_hour = 0
    /* 时 - 取值区间为[0,23] */
    val tm_mday = 0
    /* 一个月中的日期 - 取值区间为[1,31] */
    val tm_mon = 0
    /* 月份（从一月开始，0代表一月） - 取值区间为[0,11] */
    val tm_year = 0
    /* 年份，其值等于实际年份减去1900 */
    val tm_wday = 0
    /* 星期 – 取值区间为[0,6]，其中0代表星期天，1代表星期一，以此类推 */
    val tm_yday = 0
    /* 从每年的1月1日开始的天数 – 取值区间为[0,365]，其中0代表1月1日，1代表1月2日，以此类推 */
    val tm_isdst = 0 /* 夏令时标识符，实行夏令时的时候，tm_isdst为正。不实行夏令时的时候，tm_isdst为0；不了解情况时，tm_isdst()为负。*/
  }

  def Julian_Date(cdate: tm): Double = {
    /* The function Julian_Date converts a standard calendar   *//* date and time to a Julian Date. The procedure Date_Time *//* performs the inverse of this function. */ var julian_date = .0
    julian_date = Julian_Date_of_Year(cdate.tm_year) + DOY(cdate.tm_year, cdate.tm_mon, cdate.tm_mday) + Fraction_of_Day(cdate.tm_hour, cdate.tm_min, cdate.tm_sec) + 5.787037e-06 /* Round up to nearest 1 sec */
    julian_date
  }

  def Delta_ET(year: Double): Double = {
    /* The function Delta_ET has been added to allow calculations on   *//* the position of the sun.  It provides the difference between UT *//* (approximately the same as UTC) and ET (now referred to as TDT).*//* This function is based on a least squares fit of data from 1950 *//* to 1991 and will need to be updated periodically. *//* Values determined using data from 1950-1991 in the 1990
     Astronomical Almanac.  See DELTA_ET.WQ1 for details. */ var delta_et = .0
    delta_et = 26.465 + 0.747622 * (year - 1950) + 1.886913 * Math.sin(twopi * (year - 1975) / 33)
    delta_et
  }

  def ThetaG(epoch: Double, deep_arg: Predict#deep_arg_t): Double = {
    /* The function ThetaG calculates the Greenwich Mean Sidereal Time *//* for an epoch specified in the format used in the NORAD two-line *//* element sets. It has now been adapted for dates beyond the year *//* 1999, as described above. The function ThetaG_JD provides the   *//* same calculation except that it is based on an input in the     *//* form of a Julian Date. *//* Reference:  The 1992 Astronomical Almanac, page B6. */ var year = .0
    var day = .0
    var UT = .0
    var jd = .0
    var TU = .0
    var GMST = .0
    var ThetaG = .0
    /* Modification to support Y2K *//* Valid 1957 through 2056     */ year = modf_Integer(epoch * 1E-3)
    day = modf_Frac(epoch * 1E-3) * 1E3
    if (year < 57) year += 2000
    else year += 1900
    UT = modf_Frac(day)
    day = modf_Integer(day)
    jd = Julian_Date_of_Year(year) + day
    TU = (jd - 2451545.0) / 36525
    GMST = 24110.54841 + TU * (8640184.812866 + TU * (0.093104 - TU * 6.2E-6))
    GMST = Modulus(GMST + secday * omega_E * UT, secday)
    ThetaG = twopi * GMST / secday
    deep_arg.ds50 = jd - 2433281.5 + UT
    ThetaG = FMod2p(6.3003880987 * deep_arg.ds50 + 1.72944494)
    ThetaG
  }

  def ThetaG_JD(jdPara: Double): Double = {
    var jd = jdPara
    var UT = .0
    var TU = .0
    var GMST = .0
    UT = Frac(jd + 0.5)
    jd = jd - UT
    TU = (jd - 2451545.0) / 36525
    GMST = 24110.54841 + TU * (8640184.812866 + TU * (0.093104 - TU * 6.2E-6))
    GMST = Modulus(GMST + secday * omega_E * UT, secday)
    twopi * GMST / secday
  }

  def Calculate_Solar_Position(time: Double, solar_vector: Predict#vector_t): Unit = {
    /* Calculates solar position vector */ var mjd = .0
    var year = .0
    var T = .0
    var M = .0
    var L = .0
    var e = .0
    var C = .0
    var O = .0
    var Lsa = .0
    var nu = .0
    var R = .0
    var eps = .0
    mjd = time - 2415020.0
    year = 1900 + mjd / 365.25
    T = (mjd + Delta_ET(year) / secday) / 36525.0
    M = Radians(Modulus(358.47583 + Modulus(35999.04975 * T, 360.0) - (0.000150 + 0.0000033 * T) * Sqr(T), 360.0))
    L = Radians(Modulus(279.69668 + Modulus(36000.76892 * T, 360.0) + 0.0003025 * Sqr(T), 360.0))
    e = 0.01675104 - (0.0000418 + 0.000000126 * T) * T
    C = Radians((1.919460 - (0.004789 + 0.000014 * T) * T) * Math.sin(M) + (0.020094 - 0.000100 * T) * Math.sin(2 * M) + 0.000293 * Math.sin(3 * M))
    O = Radians(Modulus(259.18 - 1934.142 * T, 360.0))
    Lsa = Modulus(L + C - Radians(0.00569 - 0.00479 * Math.sin(O)), twopi)
    nu = Modulus(M + C, twopi)
    R = 1.0000002 * (1.0 - Sqr(e)) / (1.0 + e * Math.cos(nu))
    eps = Radians(23.452294 - (0.0130125 + (0.00000164 - 0.000000503 * T) * T) * T + 0.00256 * Math.cos(O))
    R = AU * R
    solar_vector.x = R * Math.cos(Lsa)
    solar_vector.y = R * Math.sin(Lsa) * Math.cos(eps)
    solar_vector.z = R * Math.sin(Lsa) * Math.sin(eps)
    solar_vector.w = R
  }

  def Sat_Eclipsed(pos: vector_t, sol: vector_t, depth: Double) = {
    /* Calculates satellite's eclipse status and depth */ var ret_depth = 0
    var sd_sun = .0
    var sd_earth = .0
    var delta = .0
    val Rho = new vector_t
    val earth = new vector_t
    /* Determine partial eclipse */ sd_earth = ArcSin(xkmper / pos.w)
    Vec_Sub(sol, pos, Rho)
    sd_sun = ArcSin(sr / Rho.w)
    Scalar_Multiply(-1, pos, earth)
    delta = Angle(sol, earth)
    ret_depth = (sd_earth - sd_sun - delta).toInt
    if (sd_earth < sd_sun) (false, ret_depth)
    else if (ret_depth >= 0) (true, ret_depth)
    else (false, ret_depth)
  }

  def select_ephemeris(tle: Predict#tle_t): Unit = {
    /* Selects the apropriate ephemeris type to be used *//* for predictions according to the data in the TLE *//* It also processes values in the tle set so that  *//* they are apropriate for the sgp4/sdp4 routines   */ var ao = .0
    var xnodp = .0
    var dd1 = .0
    var dd2 = .0
    var delo = .0
    var temp = .0
    var a1 = .0
    var del1 = .0
    var r1 = .0
    /* Preprocess tle set */ tle.xnodeo *= deg2rad
    tle.omegao *= deg2rad
    tle.xmo *= deg2rad
    tle.xincl *= deg2rad
    temp = twopi / xmnpda / xmnpda
    tle.xno = tle.xno * temp * xmnpda
    tle.xndt2o *= temp
    tle.xndd6o = tle.xndd6o * temp / xmnpda
    tle.bstar /= ae
    /* Period > 225 minutes is deep space */ dd1 = xke / tle.xno
    dd2 = tothrd
    a1 = Math.pow(dd1, dd2)
    r1 = Math.cos(tle.xincl)
    dd1 = 1.0 - tle.eo * tle.eo
    temp = ck2 * 1.5f * (r1 * r1 * 3.0 - 1.0) / Math.pow(dd1, 1.5)
    del1 = temp / (a1 * a1)
    ao = a1 * (1.0 - del1 * (tothrd * .5 + del1 * (del1 * 1.654320987654321 + 1.0)))
    delo = temp / (ao * ao)
    xnodp = tle.xno / (delo + 1.0)
    /* Select a deep-space/near-earth ephemeris */ if (twopi / xnodp / xmnpda >= 0.15625) SetFlag(DEEP_SPACE_EPHEM_FLAG)
    else ClearFlag(DEEP_SPACE_EPHEM_FLAG)
  }

  val sgp4Object = new SGP4Class

  def SGP4(tsince: Double, tle: Predict#tle_t, pos: Predict#vector_t, vel: Predict#vector_t): Unit = {
    sgp4Object.SGP4(tsince, tle, pos, vel)
  }

  class SGP4Class {
    var aodp = .0
    var aycof = .0
    var c1 = .0
    var c4 = .0
    var c5 = .0
    var cosio = .0
    var d2 = .0
    var d3 = .0
    var d4 = .0
    var delmo = .0
    var omgcof = .0
    var eta = .0
    var omgdot = .0
    var sinio = .0
    var xnodp = .0
    var sinmo = .0
    var t2cof = .0
    var t3cof = .0
    var t4cof = .0
    var t5cof = .0
    var x1mth2 = .0
    var x3thm1 = .0
    var x7thm1 = .0
    var xmcof = .0
    var xmdot = .0
    var xnodcf = .0
    var xnodot = .0
    var xlcof = .0

    def SGP4(tsince: Double, tle: Predict#tle_t, pos: Predict#vector_t, vel: Predict#vector_t): Unit = {
      /* This function is used to calculate the position and velocity *//* of near-earth (period < 225 minutes) satellites. tsince is   *//* time since epoch in minutes, tle is a pointer to a tle_t     *//* structure with Keplerian orbital elements and pos and vel    *//* are vector_t structures returning ECI satellite position and *//* velocity. Use Convert_Sat_State() to convert to km and km/s. */ var cosuk = .0
      var sinuk = .0
      var rfdotk = .0
      var vx = .0
      var vy = .0
      var vz = .0
      var ux = .0
      var uy = .0
      var uz = .0
      var xmy = .0
      var xmx = .0
      var cosnok = .0
      var sinnok = .0
      var cosik = .0
      var sinik = .0
      var rdotk = .0
      var xinck = .0
      var xnodek = .0
      var uk = .0
      var rk = .0
      var cos2u = .0
      var sin2u = .0
      var u = .0
      var sinu = .0
      var cosu = .0
      var betal = .0
      var rfdot = .0
      var rdot = .0
      var r = .0
      var pl = .0
      var elsq = .0
      var esine = .0
      var ecose = .0
      var epw = .0
      var cosepw = .0
      var x1m5th = .0
      var xhdot1 = .0
      var tfour = .0
      var sinepw = .0
      var capu = .0
      var ayn = .0
      var xlt = .0
      var aynl = .0
      var xll = .0
      var axn = .0
      var xn = .0
      var beta = .0
      var xl = .0
      var e = .0
      var a = .0
      var tcube = .0
      var delm = .0
      var delomg = .0
      var templ = .0
      var tempe = .0
      var tempa = .0
      var xnode = .0
      var tsq = .0
      var xmp = .0
      var omega = .0
      var xnoddf = .0
      var omgadf = .0
      var xmdf = .0
      var a1 = .0
      var a3ovk2 = .0
      var ao = .0
      var betao = .0
      var betao2 = .0
      var c1sq = .0
      var c2 = .0
      var c3 = .0
      var coef = .0
      var coef1 = .0
      var del1 = .0
      var delo = .0
      var eeta = .0
      var eosq = .0
      var etasq = .0
      var perigee = .0
      var pinvsq = .0
      var psisq = .0
      var qoms24 = .0
      var s4 = .0
      var temp = .0
      var temp1 = .0
      var temp2 = .0
      var temp3 = .0
      var temp4 = .0
      var temp5 = .0
      var temp6 = .0
      var theta2 = .0
      var theta4 = .0
      var tsi = .0
      var i = 0
      /* Initialization */ if (isFlagClear(SGP4_INITIALIZED_FLAG)) {
        SetFlag(SGP4_INITIALIZED_FLAG)
        /* Recover original mean motion (xnodp) and   *//* semimajor axis (aodp) from input elements. */ a1 = Math.pow(xke / tle.xno, tothrd)
        cosio = Math.cos(tle.xincl)
        theta2 = cosio * cosio
        x3thm1 = 3 * theta2 - 1.0
        eosq = tle.eo * tle.eo
        betao2 = 1.0 - eosq
        betao = Math.sqrt(betao2)
        del1 = 1.5 * ck2 * x3thm1 / (a1 * a1 * betao * betao2)
        ao = a1 * (1.0 - del1 * (0.5 * tothrd + del1 * (1.0 + 134.0 / 81.0 * del1)))
        delo = 1.5 * ck2 * x3thm1 / (ao * ao * betao * betao2)
        xnodp = tle.xno / (1.0 + delo)
        aodp = ao / (1.0 - delo)
        /* For perigee less than 220 kilometers, the "simple"     *//* flag is set and the equations are truncated to linear  *//* variation in sqrt a and quadratic variation in mean    *//* anomaly.  Also, the c3 term, the delta omega term, and *//* the delta m term are dropped.                          */ if ((aodp * (1 - tle.eo) / ae) < (220 / xkmper + ae)) SetFlag(SIMPLE_FLAG)
        else ClearFlag(SIMPLE_FLAG)
        /* For perigees below 156 km, the      *//* values of s and qoms2t are altered. */ s4 = s
        qoms24 = qoms2t
        perigee = (aodp * (1 - tle.eo) - ae) * xkmper
        if (perigee < 156.0) {
          if (perigee <= 98.0) s4 = 20
          else s4 = perigee - 78.0
          qoms24 = Math.pow((120 - s4) * ae / xkmper, 4)
          s4 = s4 / xkmper + ae
        }
        pinvsq = 1 / (aodp * aodp * betao2 * betao2)
        tsi = 1 / (aodp - s4)
        eta = aodp * tle.eo * tsi
        etasq = eta * eta
        eeta = tle.eo * eta
        psisq = Math.abs(1 - etasq)
        coef = qoms24 * Math.pow(tsi, 4)
        coef1 = coef / Math.pow(psisq, 3.5)
        c2 = coef1 * xnodp * (aodp * (1 + 1.5 * etasq + eeta * (4 + etasq)) + 0.75 * ck2 * tsi / psisq * x3thm1 * (8 + 3 * etasq * (8 + etasq)))
        c1 = tle.bstar * c2
        sinio = Math.sin(tle.xincl)
        a3ovk2 = -xj3 / ck2 * Math.pow(ae, 3)
        c3 = coef * tsi * a3ovk2 * xnodp * ae * sinio / tle.eo
        x1mth2 = 1 - theta2
        c4 = 2 * xnodp * coef1 * aodp * betao2 * (eta * (2 + 0.5 * etasq) + tle.eo * (0.5 + 2 * etasq) - 2 * ck2 * tsi / (aodp * psisq) * (-3 * x3thm1 * (1 - 2 * eeta + etasq * (1.5 - 0.5 * eeta)) + 0.75 * x1mth2 * (2 * etasq - eeta * (1 + etasq)) * Math.cos(2 * tle.omegao)))
        c5 = 2 * coef1 * aodp * betao2 * (1 + 2.75 * (etasq + eeta) + eeta * etasq)
        theta4 = theta2 * theta2
        temp1 = 3 * ck2 * pinvsq * xnodp
        temp2 = temp1 * ck2 * pinvsq
        temp3 = 1.25 * ck4 * pinvsq * pinvsq * xnodp
        xmdot = xnodp + 0.5 * temp1 * betao * x3thm1 + 0.0625 * temp2 * betao * (13 - 78 * theta2 + 137 * theta4)
        x1m5th = 1 - 5 * theta2
        omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7 - 114 * theta2 + 395 * theta4) + temp3 * (3 - 36 * theta2 + 49 * theta4)
        xhdot1 = -temp1 * cosio
        xnodot = xhdot1 + (0.5 * temp2 * (4 - 19 * theta2) + 2 * temp3 * (3 - 7 * theta2)) * cosio
        omgcof = tle.bstar * c3 * Math.cos(tle.omegao)
        xmcof = -tothrd * coef * tle.bstar * ae / eeta
        xnodcf = 3.5 * betao2 * xhdot1 * c1
        t2cof = 1.5 * c1
        xlcof = 0.125 * a3ovk2 * sinio * (3 + 5 * cosio) / (1 + cosio)
        aycof = 0.25 * a3ovk2 * sinio
        delmo = Math.pow(1 + eta * Math.cos(tle.xmo), 3)
        sinmo = Math.sin(tle.xmo)
        x7thm1 = 7 * theta2 - 1
        if (isFlagClear(SIMPLE_FLAG)) {
          c1sq = c1 * c1
          d2 = 4 * aodp * tsi * c1sq
          temp = d2 * tsi * c1 / 3
          d3 = (17 * aodp + s4) * temp
          d4 = 0.5 * temp * aodp * tsi * (221 * aodp + 31 * s4) * c1
          t3cof = d2 + 2 * c1sq
          t4cof = 0.25 * (3 * d3 + c1 * (12 * d2 + 10 * c1sq))
          t5cof = 0.2 * (3 * d4 + 12 * c1 * d3 + 6 * d2 * d2 + 15 * c1sq * (2 * d2 + c1sq))
        }
      }
      /* Update for secular gravity and atmospheric drag. */ xmdf = tle.xmo + xmdot * tsince
      omgadf = tle.omegao + omgdot * tsince
      xnoddf = tle.xnodeo + xnodot * tsince
      omega = omgadf
      xmp = xmdf
      tsq = tsince * tsince
      xnode = xnoddf + xnodcf * tsq
      tempa = 1 - c1 * tsince
      tempe = tle.bstar * c4 * tsince
      templ = t2cof * tsq
      if (isFlagClear(SIMPLE_FLAG)) {
        delomg = omgcof * tsince
        delm = xmcof * (Math.pow(1 + eta * Math.cos(xmdf), 3) - delmo)
        temp = delomg + delm
        xmp = xmdf + temp
        omega = omgadf - temp
        tcube = tsq * tsince
        tfour = tsince * tcube
        tempa = tempa - d2 * tsq - d3 * tcube - d4 * tfour
        tempe = tempe + tle.bstar * c5 * (Math.sin(xmp) - sinmo)
        templ = templ + t3cof * tcube + tfour * (t4cof + tsince * t5cof)
      }
      a = aodp * Math.pow(tempa, 2)
      e = tle.eo - tempe
      xl = xmp + omega + xnode + xnodp * templ
      beta = Math.sqrt(1 - e * e)
      xn = xke / Math.pow(a, 1.5)
      /* Long period periodics */ axn = e * Math.cos(omega)
      temp = 1 / (a * beta * beta)
      xll = temp * xlcof * axn
      aynl = temp * aycof
      xlt = xl + xll
      ayn = e * Math.sin(omega) + aynl
      /* Solve Kepler's Equation */ capu = FMod2p(xlt - xnode)
      temp2 = capu
      i = 0
      var needBreak = false
      do {
        sinepw = Math.sin(temp2)
        cosepw = Math.cos(temp2)
        temp3 = axn * sinepw
        temp4 = ayn * cosepw
        temp5 = axn * cosepw
        temp6 = ayn * sinepw
        epw = (capu - temp4 + temp3 - temp2) / (1 - temp5 - temp6) + temp2
        if (Math.abs(epw - temp2) <= e6a) needBreak = true //todo: break is not supported
        if (!needBreak) temp2 = epw
      } while ( {
        {
          i += 1;
          i - 1
        } < 10 && !needBreak
      })
      /* Short period preliminary quantities */ ecose = temp5 + temp6
      esine = temp3 - temp4
      elsq = axn * axn + ayn * ayn
      temp = 1 - elsq
      pl = a * temp
      r = a * (1 - ecose)
      temp1 = 1 / r
      rdot = xke * Math.sqrt(a) * esine * temp1
      rfdot = xke * Math.sqrt(pl) * temp1
      temp2 = a * temp1
      betal = Math.sqrt(temp)
      temp3 = 1 / (1 + betal)
      cosu = temp2 * (cosepw - axn + ayn * esine * temp3)
      sinu = temp2 * (sinepw - ayn - axn * esine * temp3)
      u = AcTan(sinu, cosu)
      sin2u = 2 * sinu * cosu
      cos2u = 2 * cosu * cosu - 1
      temp = 1 / pl
      temp1 = ck2 * temp
      temp2 = temp1 * temp
      /* Update for short periodics */ rk = r * (1 - 1.5 * temp2 * betal * x3thm1) + 0.5 * temp1 * x1mth2 * cos2u
      uk = u - 0.25 * temp2 * x7thm1 * sin2u
      xnodek = xnode + 1.5 * temp2 * cosio * sin2u
      xinck = tle.xincl + 1.5 * temp2 * cosio * sinio * cos2u
      rdotk = rdot - xn * temp1 * x1mth2 * sin2u
      rfdotk = rfdot + xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1)
      /* Orientation vectors */ sinuk = Math.sin(uk)
      cosuk = Math.cos(uk)
      sinik = Math.sin(xinck)
      cosik = Math.cos(xinck)
      sinnok = Math.sin(xnodek)
      cosnok = Math.cos(xnodek)
      xmx = -sinnok * cosik
      xmy = cosnok * cosik
      ux = xmx * sinuk + cosnok * cosuk
      uy = xmy * sinuk + sinnok * cosuk
      uz = sinik * sinuk
      vx = xmx * cosuk - cosnok * sinuk
      vy = xmy * cosuk - sinnok * sinuk
      vz = sinik * cosuk
      /* Position and velocity */ pos.x = rk * ux
      pos.y = rk * uy
      pos.z = rk * uz
      vel.x = rdotk * ux + rfdotk * vx
      vel.y = rdotk * uy + rfdotk * vy
      vel.z = rdotk * uz + rfdotk * vz
      /* Phase in radians */ phase = xlt - xnode - omgadf + twopi
      if (phase < 0.0) phase += twopi
      phase = FMod2p(phase)
    }
  }

  val deepObject = new DeepClass

  def Deep(ientry: Int, tle: tle_t, deep_arg: deep_arg_t): Unit = {
    deepObject.Deep(ientry, tle, deep_arg)
  }

  class DeepClass {
    var thgr = .0
    var xnq = .0
    var xqncl = .0
    var omegaq = .0
    var zmol = .0
    var zmos = .0
    var savtsn = .0
    var ee2 = .0
    var e3 = .0
    var xi2 = .0
    var xl2 = .0
    var xl3 = .0
    var xl4 = .0
    var xgh2 = .0
    var xgh3 = .0
    var xgh4 = .0
    var xh2 = .0
    var xh3 = .0
    var sse = .0
    var ssi = .0
    var ssg = .0
    var xi3 = .0
    var se2 = .0
    var si2 = .0
    var sl2 = .0
    var sgh2 = .0
    var sh2 = .0
    var se3 = .0
    var si3 = .0
    var sl3 = .0
    var sgh3 = .0
    var sh3 = .0
    var sl4 = .0
    var sgh4 = .0
    var ssl = .0
    var ssh = .0
    var d3210 = .0
    var d3222 = .0
    var d4410 = .0
    var d4422 = .0
    var d5220 = .0
    var d5232 = .0
    var d5421 = .0
    var d5433 = .0
    var del1 = .0
    var del2 = .0
    var del3 = .0
    var fasx2 = .0
    var fasx4 = .0
    var fasx6 = .0
    var xlamo = .0
    var xfact = .0
    var xni = .0
    var atime = .0
    var stepp = .0
    var stepn = .0
    var step2 = .0
    var preep = .0
    var pl = .0
    var sghs = .0
    var xli = .0
    var d2201 = .0
    var d2211 = .0
    var sghl = .0
    var sh1 = .0
    var pinc = .0
    var pe = .0
    var shs = .0
    var zsingl = .0
    var zcosgl = .0
    var zsinhl = .0
    var zcoshl = .0
    var zsinil = .0
    var zcosil = .0

    def Deep(ientry: Int, tle: tle_t, deep_arg: deep_arg_t): Unit = {
      /* This function is used by SDP4 to add lunar and solar *//* perturbation effects to deep-space orbit objects.    */ var a1 = .0
      var a2 = .0
      var a3 = .0
      var a4 = .0
      var a5 = .0
      var a6 = .0
      var a7 = .0
      var a8 = .0
      var a9 = .0
      var a10 = .0
      var ainv2 = .0
      var alfdp = .0
      var aqnv = .0
      var sgh = .0
      var sini2 = .0
      var sinis = .0
      var sinok = .0
      var sh = .0
      var si = .0
      var sil = .0
      var day = .0
      var betdp = .0
      var dalf = .0
      var bfact = .0
      var c = .0
      var cc = .0
      var cosis = .0
      var cosok = .0
      var cosq = .0
      var ctem = .0
      var f322 = .0
      var zx = .0
      var zy = .0
      var dbet = .0
      var dls = .0
      var eoc = .0
      var eq = .0
      var f2 = .0
      var f220 = .0
      var f221 = .0
      var f3 = .0
      var f311 = .0
      var f321 = .0
      var xnoh = .0
      var f330 = .0
      var f441 = .0
      var f442 = .0
      var f522 = .0
      var f523 = .0
      var f542 = .0
      var f543 = .0
      var g200 = .0
      var g201 = .0
      var g211 = .0
      var pgh = .0
      var ph = .0
      var s1 = .0
      var s2 = .0
      var s3 = .0
      var s4 = .0
      var s5 = .0
      var s6 = .0
      var s7 = .0
      var se = .0
      var sel = .0
      var ses = .0
      var xls = .0
      var g300 = .0
      var g310 = .0
      var g322 = .0
      var g410 = .0
      var g422 = .0
      var g520 = .0
      var g521 = .0
      var g532 = .0
      var g533 = .0
      var gam = .0
      var sinq = .0
      var sinzf = .0
      var sis = .0
      var sl = .0
      var sll = .0
      var sls = .0
      var stem = .0
      var temp = .0
      var temp1 = .0
      var x1 = .0
      var x2 = .0
      var x2li = .0
      var x2omi = .0
      var x3 = .0
      var x4 = .0
      var x5 = .0
      var x6 = .0
      var x7 = .0
      var x8 = .0
      var xl = .0
      var xldot = .0
      var xmao = .0
      var xnddt = .0
      var xndot = .0
      var xno2 = .0
      var xnodce = .0
      var xnoi = .0
      var xomi = .0
      var xpidot = .0
      var z1 = .0
      var z11 = .0
      var z12 = .0
      var z13 = .0
      var z2 = .0
      var z21 = .0
      var z22 = .0
      var z23 = .0
      var z3 = .0
      var z31 = .0
      var z32 = .0
      var z33 = .0
      var ze = .0
      var zf = .0
      var zm = .0
      var zn = .0
      var zsing = .0
      var zsinh = .0
      var zsini = .0
      var zcosg = .0
      var zcosh = .0
      var zcosi = .0
      var delt = 0
      var ft = 0
      if (ientry == dpinit) {
        /* Entrance for deep space initialization */ thgr = ThetaG(tle.epoch, deep_arg)
        eq = tle.eo
        xnq = deep_arg.xnodp
        aqnv = 1 / deep_arg.aodp
        xqncl = tle.xincl
        xmao = tle.xmo
        xpidot = deep_arg.omgdot + deep_arg.xnodot
        sinq = Math.sin(tle.xnodeo)
        cosq = Math.cos(tle.xnodeo)
        omegaq = tle.omegao
        /* Initialize lunar solar terms */ day = deep_arg.ds50 + 18261.5
        /* Days since 1900 Jan 0.5 */ if (day != preep) {
          preep = day
          xnodce = 4.5236020 - 9.2422029E-4 * day
          stem = Math.sin(xnodce)
          ctem = Math.cos(xnodce)
          zcosil = 0.91375164 - 0.03568096 * ctem
          zsinil = Math.sqrt(1 - zcosil * zcosil)
          zsinhl = 0.089683511 * stem / zsinil
          zcoshl = Math.sqrt(1 - zsinhl * zsinhl)
          c = 4.7199672 + 0.22997150 * day
          gam = 5.8351514 + 0.0019443680 * day
          zmol = FMod2p(c - gam)
          zx = 0.39785416 * stem / zsinil
          zy = zcoshl * ctem + 0.91744867 * zsinhl * stem
          zx = AcTan(zx, zy)
          zx = gam + zx - xnodce
          zcosgl = Math.cos(zx)
          zsingl = Math.sin(zx)
          zmos = 6.2565837 + 0.017201977 * day
          zmos = FMod2p(zmos)
        }
        /* Do solar terms */ savtsn = 1E20
        zcosg = zcosgs
        zsing = zsings
        zcosi = zcosis
        zsini = zsinis
        zcosh = cosq
        zsinh = sinq
        cc = c1ss
        zn = zns
        ze = zes
        xnoi = 1 / xnq
        /* Loop breaks when Solar terms are done a second *//* time, after Lunar terms are initialized        */

        var needBreak = false
        while (!needBreak) {
          /* Solar terms done again after Lunar terms are done */ a1 = zcosg * zcosh + zsing * zcosi * zsinh
          a3 = -zsing * zcosh + zcosg * zcosi * zsinh
          a7 = -zcosg * zsinh + zsing * zcosi * zcosh
          a8 = zsing * zsini
          a9 = zsing * zsinh + zcosg * zcosi * zcosh
          a10 = zcosg * zsini
          a2 = deep_arg.cosio * a7 + deep_arg.sinio * a8
          a4 = deep_arg.cosio * a9 + deep_arg.sinio * a10
          a5 = -deep_arg.sinio * a7 + deep_arg.cosio * a8
          a6 = -deep_arg.sinio * a9 + deep_arg.cosio * a10
          x1 = a1 * deep_arg.cosg + a2 * deep_arg.sing
          x2 = a3 * deep_arg.cosg + a4 * deep_arg.sing
          x3 = -a1 * deep_arg.sing + a2 * deep_arg.cosg
          x4 = -a3 * deep_arg.sing + a4 * deep_arg.cosg
          x5 = a5 * deep_arg.sing
          x6 = a6 * deep_arg.sing
          x7 = a5 * deep_arg.cosg
          x8 = a6 * deep_arg.cosg
          z31 = 12 * x1 * x1 - 3 * x3 * x3
          z32 = 24 * x1 * x2 - 6 * x3 * x4
          z33 = 12 * x2 * x2 - 3 * x4 * x4
          z1 = 3 * (a1 * a1 + a2 * a2) + z31 * deep_arg.eosq
          z2 = 6 * (a1 * a3 + a2 * a4) + z32 * deep_arg.eosq
          z3 = 3 * (a3 * a3 + a4 * a4) + z33 * deep_arg.eosq
          z11 = -6 * a1 * a5 + deep_arg.eosq * (-24 * x1 * x7 - 6 * x3 * x5)
          z12 = -6 * (a1 * a6 + a3 * a5) + deep_arg.eosq * (-24 * (x2 * x7 + x1 * x8) - 6 * (x3 * x6 + x4 * x5))
          z13 = -6 * a3 * a6 + deep_arg.eosq * (-24 * x2 * x8 - 6 * x4 * x6)
          z21 = 6 * a2 * a5 + deep_arg.eosq * (24 * x1 * x5 - 6 * x3 * x7)
          z22 = 6 * (a4 * a5 + a2 * a6) + deep_arg.eosq * (24 * (x2 * x5 + x1 * x6) - 6 * (x4 * x7 + x3 * x8))
          z23 = 6 * a4 * a6 + deep_arg.eosq * (24 * x2 * x6 - 6 * x4 * x8)
          z1 = z1 + z1 + deep_arg.betao2 * z31
          z2 = z2 + z2 + deep_arg.betao2 * z32
          z3 = z3 + z3 + deep_arg.betao2 * z33
          s3 = cc * xnoi
          s2 = -0.5 * s3 / deep_arg.betao
          s4 = s3 * deep_arg.betao
          s1 = -15 * eq * s4
          s5 = x1 * x3 + x2 * x4
          s6 = x2 * x3 + x1 * x4
          s7 = x2 * x4 - x1 * x3
          se = s1 * zn * s5
          si = s2 * zn * (z11 + z13)
          sl = -zn * s3 * (z1 + z3 - 14 - 6 * deep_arg.eosq)
          sgh = s4 * zn * (z31 + z33 - 6)
          sh = -zn * s2 * (z21 + z23)
          if (xqncl < 5.2359877E-2) sh = 0
          ee2 = 2 * s1 * s6
          e3 = 2 * s1 * s7
          xi2 = 2 * s2 * z12
          xi3 = 2 * s2 * (z13 - z11)
          xl2 = -2 * s3 * z2
          xl3 = -2 * s3 * (z3 - z1)
          xl4 = -2 * s3 * (-21 - 9 * deep_arg.eosq) * ze
          xgh2 = 2 * s4 * z32
          xgh3 = 2 * s4 * (z33 - z31)
          xgh4 = -18 * s4 * ze
          xh2 = -2 * s2 * z22
          xh3 = -2 * s2 * (z23 - z21)
          if (isFlagSet(LUNAR_TERMS_DONE_FLAG)) needBreak = true //todo: break is not supported
          if (!needBreak) {
            /* Do lunar terms */ sse = se
            ssi = si
            ssl = sl
            ssh = sh / deep_arg.sinio
            ssg = sgh - deep_arg.cosio * ssh
            se2 = ee2
            si2 = xi2
            sl2 = xl2
            sgh2 = xgh2
            sh2 = xh2
            se3 = e3
            si3 = xi3
            sl3 = xl3
            sgh3 = xgh3
            sh3 = xh3
            sl4 = xl4
            sgh4 = xgh4
            zcosg = zcosgl
            zsing = zsingl
            zcosi = zcosil
            zsini = zsinil
            zcosh = zcoshl * cosq + zsinhl * sinq
            zsinh = sinq * zcoshl - cosq * zsinhl
            zn = znl
            cc = c1l
            ze = zel
            SetFlag(LUNAR_TERMS_DONE_FLAG)
          }
        }
        sse = sse + se
        ssi = ssi + si
        ssl = ssl + sl
        ssg = ssg + sgh - deep_arg.cosio / deep_arg.sinio * sh
        ssh = ssh + sh / deep_arg.sinio
        /* Geopotential resonance initialization for 12 hour orbits */ ClearFlag(RESONANCE_FLAG)
        ClearFlag(SYNCHRONOUS_FLAG)
        if (!((xnq < 0.0052359877) && (xnq > 0.0034906585))) {
          if ((xnq < 0.00826) || (xnq > 0.00924)) return
          if (eq < 0.5) return
          SetFlag(RESONANCE_FLAG)
          eoc = eq * deep_arg.eosq
          g201 = -0.306 - (eq - 0.64) * 0.440
          if (eq <= 0.65) {
            g211 = 3.616 - 13.247 * eq + 16.290 * deep_arg.eosq
            g310 = -19.302 + 117.390 * eq - 228.419 * deep_arg.eosq + 156.591 * eoc
            g322 = -18.9068 + 109.7927 * eq - 214.6334 * deep_arg.eosq + 146.5816 * eoc
            g410 = -41.122 + 242.694 * eq - 471.094 * deep_arg.eosq + 313.953 * eoc
            g422 = -146.407 + 841.880 * eq - 1629.014 * deep_arg.eosq + 1083.435 * eoc
            g520 = -532.114 + 3017.977 * eq - 5740 * deep_arg.eosq + 3708.276 * eoc
          }
          else {
            g211 = -72.099 + 331.819 * eq - 508.738 * deep_arg.eosq + 266.724 * eoc
            g310 = -346.844 + 1582.851 * eq - 2415.925 * deep_arg.eosq + 1246.113 * eoc
            g322 = -342.585 + 1554.908 * eq - 2366.899 * deep_arg.eosq + 1215.972 * eoc
            g410 = -1052.797 + 4758.686 * eq - 7193.992 * deep_arg.eosq + 3651.957 * eoc
            g422 = -3581.69 + 16178.11 * eq - 24462.77 * deep_arg.eosq + 12422.52 * eoc
            if (eq <= 0.715) g520 = 1464.74 - 4664.75 * eq + 3763.64 * deep_arg.eosq
            else g520 = -5149.66 + 29936.92 * eq - 54087.36 * deep_arg.eosq + 31324.56 * eoc
          }
          if (eq < 0.7) {
            g533 = -919.2277 + 4988.61 * eq - 9064.77 * deep_arg.eosq + 5542.21 * eoc
            g521 = -822.71072 + 4568.6173 * eq - 8491.4146 * deep_arg.eosq + 5337.524 * eoc
            g532 = -853.666 + 4690.25 * eq - 8624.77 * deep_arg.eosq + 5341.4 * eoc
          }
          else {
            g533 = -37995.78 + 161616.52 * eq - 229838.2 * deep_arg.eosq + 109377.94 * eoc
            g521 = -51752.104 + 218913.95 * eq - 309468.16 * deep_arg.eosq + 146349.42 * eoc
            g532 = -40023.88 + 170470.89 * eq - 242699.48 * deep_arg.eosq + 115605.82 * eoc
          }
          sini2 = deep_arg.sinio * deep_arg.sinio
          f220 = 0.75 * (1 + 2 * deep_arg.cosio + deep_arg.theta2)
          f221 = 1.5 * sini2
          f321 = 1.875 * deep_arg.sinio * (1 - 2 * deep_arg.cosio - 3 * deep_arg.theta2)
          f322 = -1.875 * deep_arg.sinio * (1 + 2 * deep_arg.cosio - 3 * deep_arg.theta2)
          f441 = 35 * sini2 * f220
          f442 = 39.3750 * sini2 * sini2
          f522 = 9.84375 * deep_arg.sinio * (sini2 * (1 - 2 * deep_arg.cosio - 5 * deep_arg.theta2) + 0.33333333 * (-2 + 4 * deep_arg.cosio + 6 * deep_arg.theta2))
          f523 = deep_arg.sinio * (4.92187512 * sini2 * (-2 - 4 * deep_arg.cosio + 10 * deep_arg.theta2) + 6.56250012 * (1 + 2 * deep_arg.cosio - 3 * deep_arg.theta2))
          f542 = 29.53125 * deep_arg.sinio * (2 - 8 * deep_arg.cosio + deep_arg.theta2 * (-12 + 8 * deep_arg.cosio + 10 * deep_arg.theta2))
          f543 = 29.53125 * deep_arg.sinio * (-2 - 8 * deep_arg.cosio + deep_arg.theta2 * (12 + 8 * deep_arg.cosio - 10 * deep_arg.theta2))
          xno2 = xnq * xnq
          ainv2 = aqnv * aqnv
          temp1 = 3 * xno2 * ainv2
          temp = temp1 * root22
          d2201 = temp * f220 * g201
          d2211 = temp * f221 * g211
          temp1 = temp1 * aqnv
          temp = temp1 * root32
          d3210 = temp * f321 * g310
          d3222 = temp * f322 * g322
          temp1 = temp1 * aqnv
          temp = 2 * temp1 * root44
          d4410 = temp * f441 * g410
          d4422 = temp * f442 * g422
          temp1 = temp1 * aqnv
          temp = temp1 * root52
          d5220 = temp * f522 * g520
          d5232 = temp * f523 * g532
          temp = 2 * temp1 * root54
          d5421 = temp * f542 * g521
          d5433 = temp * f543 * g533
          xlamo = xmao + tle.xnodeo + tle.xnodeo - thgr - thgr
          bfact = deep_arg.xmdot + deep_arg.xnodot + deep_arg.xnodot - thdt - thdt
          bfact = bfact + ssl + ssh + ssh
        }
        else {
          SetFlag(RESONANCE_FLAG)
          SetFlag(SYNCHRONOUS_FLAG)
          /* Synchronous resonance terms initialization */ g200 = 1 + deep_arg.eosq * (-2.5 + 0.8125 * deep_arg.eosq)
          g310 = 1 + 2 * deep_arg.eosq
          g300 = 1 + deep_arg.eosq * (-6 + 6.60937 * deep_arg.eosq)
          f220 = 0.75 * (1 + deep_arg.cosio) * (1 + deep_arg.cosio)
          f311 = 0.9375 * deep_arg.sinio * deep_arg.sinio * (1 + 3 * deep_arg.cosio) - 0.75 * (1 + deep_arg.cosio)
          f330 = 1 + deep_arg.cosio
          f330 = 1.875 * f330 * f330 * f330
          del1 = 3 * xnq * xnq * aqnv * aqnv
          del2 = 2 * del1 * f220 * g200 * q22
          del3 = 3 * del1 * f330 * g300 * q33 * aqnv
          del1 = del1 * f311 * g310 * q31 * aqnv
          fasx2 = 0.13130908
          fasx4 = 2.8843198
          fasx6 = 0.37448087
          xlamo = xmao + tle.xnodeo + tle.omegao - thgr
          bfact = deep_arg.xmdot + xpidot - thdt
          bfact = bfact + ssl + ssg + ssh
        }
        xfact = bfact - xnq
        /* Initialize integrator */ xli = xlamo
        xni = xnq
        atime = 0
        stepp = 720
        stepn = -720
        step2 = 259200
      }
      if (ientry == dpsec) {
        /* Entrance for deep space secular effects */ deep_arg.xll = deep_arg.xll + ssl * deep_arg.t
        deep_arg.omgadf = deep_arg.omgadf + ssg * deep_arg.t
        deep_arg.xnode = deep_arg.xnode + ssh * deep_arg.t
        deep_arg.em = tle.eo + sse * deep_arg.t
        deep_arg.xinc = tle.xincl + ssi * deep_arg.t
        if (deep_arg.xinc < 0) {
          deep_arg.xinc = -deep_arg.xinc
          deep_arg.xnode = deep_arg.xnode + Math.PI
          deep_arg.omgadf = deep_arg.omgadf - Math.PI
        }
        if (isFlagClear(RESONANCE_FLAG)) return
        do {
          if ((atime == 0) || ((deep_arg.t >= 0) && (atime < 0)) || ((deep_arg.t < 0) && (atime >= 0))) {
            /* Epoch restart */ if (deep_arg.t >= 0) delt = stepp.toInt
            else delt = stepn.toInt
            atime = 0
            xni = xnq
            xli = xlamo
          }
          else if (Math.abs(deep_arg.t) >= Math.abs(atime)) if (deep_arg.t > 0) delt = stepp.toInt
          else delt = stepn.toInt
          do {
            if (Math.abs(deep_arg.t - atime) >= stepp) {
              SetFlag(DO_LOOP_FLAG)
              ClearFlag(EPOCH_RESTART_FLAG)
            }
            else {
              ft = (deep_arg.t - atime).toInt
              ClearFlag(DO_LOOP_FLAG)
            }
            if (Math.abs(deep_arg.t) < Math.abs(atime)) {
              if (deep_arg.t >= 0) delt = stepn.toInt
              else delt = stepp.toInt
              SetFlag(DO_LOOP_FLAG | EPOCH_RESTART_FLAG)
            }
            /* Dot terms calculated */ if (isFlagSet(SYNCHRONOUS_FLAG)) {
              xndot = del1 * Math.sin(xli - fasx2) + del2 * Math.sin(2 * (xli - fasx4)) + del3 * Math.sin(3 * (xli - fasx6))
              xnddt = del1 * Math.cos(xli - fasx2) + 2 * del2 * Math.cos(2 * (xli - fasx4)) + 3 * del3 * Math.cos(3 * (xli - fasx6))
            }
            else {
              xomi = omegaq + deep_arg.omgdot * atime
              x2omi = xomi + xomi
              x2li = xli + xli
              xndot = d2201 * Math.sin(x2omi + xli - g22) + d2211 * Math.sin(xli - g22) + d3210 * Math.sin(xomi + xli - g32) + d3222 * Math.sin(-xomi + xli - g32) + d4410 * Math.sin(x2omi + x2li - g44) + d4422 * Math.sin(x2li - g44) + d5220 * Math.sin(xomi + xli - g52) + d5232 * Math.sin(-xomi + xli - g52) + d5421 * Math.sin(xomi + x2li - g54) + d5433 * Math.sin(-xomi + x2li - g54)
              xnddt = d2201 * Math.cos(x2omi + xli - g22) + d2211 * Math.cos(xli - g22) + d3210 * Math.cos(xomi + xli - g32) + d3222 * Math.cos(-xomi + xli - g32) + d5220 * Math.cos(xomi + xli - g52) + d5232 * Math.cos(-xomi + xli - g52) + 2 * (d4410 * Math.cos(x2omi + x2li - g44) + d4422 * Math.cos(x2li - g44) + d5421 * Math.cos(xomi + x2li - g54) + d5433 * Math.cos(-xomi + x2li - g54))
            }
            xldot = xni + xfact
            xnddt = xnddt * xldot
            if (isFlagSet(DO_LOOP_FLAG)) {
              xli = xli + xldot * delt + xndot * step2
              xni = xni + xndot * delt + xnddt * step2
              atime = atime + delt
            }
          } while ( {
            isFlagSet(DO_LOOP_FLAG) && isFlagClear(EPOCH_RESTART_FLAG)
          })
        } while ( {
          isFlagSet(DO_LOOP_FLAG) && isFlagSet(EPOCH_RESTART_FLAG)
        })
        deep_arg.xn = xni + xndot * ft + xnddt * ft * ft * 0.5
        xl = xli + xldot * ft + xndot * ft * ft * 0.5
        temp = -deep_arg.xnode + thgr + deep_arg.t * thdt
        if (isFlagClear(SYNCHRONOUS_FLAG)) deep_arg.xll = xl + temp + temp
        else deep_arg.xll = xl - deep_arg.omgadf + temp
      }
      if (ientry == dpper) {
        /* Entrance for lunar-solar periodics */ sinis = Math.sin(deep_arg.xinc)
        cosis = Math.cos(deep_arg.xinc)
        if (Math.abs(savtsn - deep_arg.t) >= 30) {
          savtsn = deep_arg.t
          zm = zmos + zns * deep_arg.t
          zf = zm + 2 * zes * Math.sin(zm)
          sinzf = Math.sin(zf)
          f2 = 0.5 * sinzf * sinzf - 0.25
          f3 = -0.5 * sinzf * Math.cos(zf)
          ses = se2 * f2 + se3 * f3
          sis = si2 * f2 + si3 * f3
          sls = sl2 * f2 + sl3 * f3 + sl4 * sinzf
          sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf
          shs = sh2 * f2 + sh3 * f3
          zm = zmol + znl * deep_arg.t
          zf = zm + 2 * zel * Math.sin(zm)
          sinzf = Math.sin(zf)
          f2 = 0.5 * sinzf * sinzf - 0.25
          f3 = -0.5 * sinzf * Math.cos(zf)
          sel = ee2 * f2 + e3 * f3
          sil = xi2 * f2 + xi3 * f3
          sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf
          sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf
          sh1 = xh2 * f2 + xh3 * f3
          pe = ses + sel
          pinc = sis + sil
          pl = sls + sll
        }
        pgh = sghs + sghl
        ph = shs + sh1
        deep_arg.xinc = deep_arg.xinc + pinc
        deep_arg.em = deep_arg.em + pe
        if (xqncl >= 0.2) {
          /* Apply periodics directly */ ph = ph / deep_arg.sinio
          pgh = pgh - deep_arg.cosio * ph
          deep_arg.omgadf = deep_arg.omgadf + pgh
          deep_arg.xnode = deep_arg.xnode + ph
          deep_arg.xll = deep_arg.xll + pl
        }
        else {
          /* Apply periodics with Lyddane modification */ sinok = Math.sin(deep_arg.xnode)
          cosok = Math.cos(deep_arg.xnode)
          alfdp = sinis * sinok
          betdp = sinis * cosok
          dalf = ph * cosok + pinc * cosis * sinok
          dbet = -ph * sinok + pinc * cosis * cosok
          alfdp = alfdp + dalf
          betdp = betdp + dbet
          deep_arg.xnode = FMod2p(deep_arg.xnode)
          xls = deep_arg.xll + deep_arg.omgadf + cosis * deep_arg.xnode
          dls = pl + pgh - pinc * deep_arg.xnode * sinis
          xls = xls + dls
          xnoh = deep_arg.xnode
          deep_arg.xnode = AcTan(alfdp, betdp)
          /* This is a patch to Lyddane modification *//* suggested by Rob Matson. */ if (Math.abs(xnoh - deep_arg.xnode) > Math.PI) if (deep_arg.xnode < xnoh) deep_arg.xnode += twopi
          else deep_arg.xnode -= twopi
          deep_arg.xll = deep_arg.xll + pl
          deep_arg.omgadf = xls - deep_arg.xll - Math.cos(deep_arg.xinc) * deep_arg.xnode
        }
      }
    }
  }

  val sdp4Object = new SDP4Class

  def SDP4(tsince: Double, tle: tle_t, pos: vector_t, vel: vector_t): Unit = {
    sdp4Object.SDP4(tsince, tle, pos, vel)
  }

  class SDP4Class {
    /* This function is used to calculate the position and velocity *//* of deep-space (period > 225 minutes) satellites. tsince is   *//* time since epoch in minutes, tle is a pointer to a tle_t     *//* structure with Keplerian orbital elements and pos and vel    *//* are vector_t structures returning ECI satellite position and *//* velocity. Use Convert_Sat_State() to convert to km and km/s. */ var x3thm1 = .0
    var c1 = .0
    var x1mth2 = .0
    var c4 = .0
    var xnodcf = .0
    var t2cof = .0
    var xlcof = .0
    var aycof = .0
    var x7thm1 = .0
    val deep_arg: deep_arg_t = null

    def SDP4(tsince: Double, tle: tle_t, pos: vector_t, vel: vector_t): Unit = {
      var i = 0
      var a = .0
      var axn = .0
      var ayn = .0
      var aynl = .0
      var beta = .0
      var betal = .0
      var capu = .0
      var cos2u = .0
      var cosepw = .0
      var cosik = .0
      var cosnok = .0
      var cosu = .0
      var cosuk = .0
      var ecose = .0
      var elsq = .0
      var epw = .0
      var esine = .0
      var pl = .0
      var theta4 = .0
      var rdot = .0
      var rdotk = .0
      var rfdot = .0
      var rfdotk = .0
      var rk = .0
      var sin2u = .0
      var sinepw = .0
      var sinik = .0
      var sinnok = .0
      var sinu = .0
      var sinuk = .0
      var tempe = .0
      var templ = .0
      var tsq = .0
      var u = .0
      var uk = .0
      var ux = .0
      var uy = .0
      var uz = .0
      var vx = .0
      var vy = .0
      var vz = .0
      var xinck = .0
      var xl = .0
      var xlt = .0
      var xmam = .0
      var xmdf = .0
      var xmx = .0
      var xmy = .0
      var xnoddf = .0
      var xnodek = .0
      var xll = .0
      var a1 = .0
      var a3ovk2 = .0
      var ao = .0
      var c2 = .0
      var coef = .0
      var coef1 = .0
      var x1m5th = .0
      var xhdot1 = .0
      var del1 = .0
      var r = .0
      var delo = .0
      var eeta = .0
      var eta = .0
      var etasq = .0
      var perigee = .0
      var psisq = .0
      var tsi = .0
      var qoms24 = .0
      var s4 = .0
      var pinvsq = .0
      var temp = .0
      var tempa = .0
      var temp1 = .0
      var temp2 = .0
      var temp3 = .0
      var temp4 = .0
      var temp5 = .0
      var temp6 = .0
      val bx = .0
      val by = .0
      val bz = .0
      val cx = .0
      val cy = .0
      val cz = .0
      /* Initialization */ if (isFlagClear(SDP4_INITIALIZED_FLAG)) {
        SetFlag(SDP4_INITIALIZED_FLAG)
        /* Recover original mean motion (xnodp) and   *//* semimajor axis (aodp) from input elements. */ a1 = Math.pow(xke / tle.xno, tothrd)
        deep_arg.cosio = Math.cos(tle.xincl)
        deep_arg.theta2 = deep_arg.cosio * deep_arg.cosio
        x3thm1 = 3 * deep_arg.theta2 - 1
        deep_arg.eosq = tle.eo * tle.eo
        deep_arg.betao2 = 1 - deep_arg.eosq
        deep_arg.betao = Math.sqrt(deep_arg.betao2)
        del1 = 1.5 * ck2 * x3thm1 / (a1 * a1 * deep_arg.betao * deep_arg.betao2)
        ao = a1 * (1 - del1 * (0.5 * tothrd + del1 * (1 + 134 / 81 * del1)))
        delo = 1.5 * ck2 * x3thm1 / (ao * ao * deep_arg.betao * deep_arg.betao2)
        deep_arg.xnodp = tle.xno / (1 + delo)
        deep_arg.aodp = ao / (1 - delo)
        /* For perigee below 156 km, the values *//* of s and qoms2t are altered.         */ s4 = s
        qoms24 = qoms2t
        perigee = (deep_arg.aodp * (1 - tle.eo) - ae) * xkmper
        if (perigee < 156.0) {
          if (perigee <= 98.0) s4 = 20.0
          else s4 = perigee - 78.0
          qoms24 = Math.pow((120 - s4) * ae / xkmper, 4)
          s4 = s4 / xkmper + ae
        }
        pinvsq = 1 / (deep_arg.aodp * deep_arg.aodp * deep_arg.betao2 * deep_arg.betao2)
        deep_arg.sing = Math.sin(tle.omegao)
        deep_arg.cosg = Math.cos(tle.omegao)
        tsi = 1 / (deep_arg.aodp - s4)
        eta = deep_arg.aodp * tle.eo * tsi
        etasq = eta * eta
        eeta = tle.eo * eta
        psisq = Math.abs(1 - etasq)
        coef = qoms24 * Math.pow(tsi, 4)
        coef1 = coef / Math.pow(psisq, 3.5)
        c2 = coef1 * deep_arg.xnodp * (deep_arg.aodp * (1 + 1.5 * etasq + eeta * (4 + etasq)) + 0.75 * ck2 * tsi / psisq * x3thm1 * (8 + 3 * etasq * (8 + etasq)))
        c1 = tle.bstar * c2
        deep_arg.sinio = Math.sin(tle.xincl)
        a3ovk2 = -xj3 / ck2 * Math.pow(ae, 3)
        x1mth2 = 1 - deep_arg.theta2
        c4 = 2 * deep_arg.xnodp * coef1 * deep_arg.aodp * deep_arg.betao2 * (eta * (2 + 0.5 * etasq) + tle.eo * (0.5 + 2 * etasq) - 2 * ck2 * tsi / (deep_arg.aodp * psisq) * (-3 * x3thm1 * (1 - 2 * eeta + etasq * (1.5 - 0.5 * eeta)) + 0.75 * x1mth2 * (2 * etasq - eeta * (1 + etasq)) * Math.cos(2 * tle.omegao)))
        theta4 = deep_arg.theta2 * deep_arg.theta2
        temp1 = 3 * ck2 * pinvsq * deep_arg.xnodp
        temp2 = temp1 * ck2 * pinvsq
        temp3 = 1.25 * ck4 * pinvsq * pinvsq * deep_arg.xnodp
        deep_arg.xmdot = deep_arg.xnodp + 0.5 * temp1 * deep_arg.betao * x3thm1 + 0.0625 * temp2 * deep_arg.betao * (13 - 78 * deep_arg.theta2 + 137 * theta4)
        x1m5th = 1 - 5 * deep_arg.theta2
        deep_arg.omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7 - 114 * deep_arg.theta2 + 395 * theta4) + temp3 * (3 - 36 * deep_arg.theta2 + 49 * theta4)
        xhdot1 = -temp1 * deep_arg.cosio
        deep_arg.xnodot = xhdot1 + (0.5 * temp2 * (4 - 19 * deep_arg.theta2) + 2 * temp3 * (3 - 7 * deep_arg.theta2)) * deep_arg.cosio
        xnodcf = 3.5 * deep_arg.betao2 * xhdot1 * c1
        t2cof = 1.5 * c1
        xlcof = 0.125 * a3ovk2 * deep_arg.sinio * (3 + 5 * deep_arg.cosio) / (1 + deep_arg.cosio)
        aycof = 0.25 * a3ovk2 * deep_arg.sinio
        x7thm1 = 7 * deep_arg.theta2 - 1
        /* initialize Deep() */ Deep(dpinit, tle, deep_arg)
      }
      /* Update for secular gravity and atmospheric drag */ xmdf = tle.xmo + deep_arg.xmdot * tsince
      deep_arg.omgadf = tle.omegao + deep_arg.omgdot * tsince
      xnoddf = tle.xnodeo + deep_arg.xnodot * tsince
      tsq = tsince * tsince
      deep_arg.xnode = xnoddf + xnodcf * tsq
      tempa = 1 - c1 * tsince
      tempe = tle.bstar * c4 * tsince
      templ = t2cof * tsq
      deep_arg.xn = deep_arg.xnodp
      /* Update for deep-space secular effects */ deep_arg.xll = xmdf
      deep_arg.t = tsince
      Deep(dpsec, tle, deep_arg)
      xmdf = deep_arg.xll
      a = Math.pow(xke / deep_arg.xn, tothrd) * tempa * tempa
      deep_arg.em = deep_arg.em - tempe
      xmam = xmdf + deep_arg.xnodp * templ
      /* Update for deep-space periodic effects */ deep_arg.xll = xmam
      Deep(dpper, tle, deep_arg)
      xmam = deep_arg.xll
      xl = xmam + deep_arg.omgadf + deep_arg.xnode
      beta = Math.sqrt(1 - deep_arg.em * deep_arg.em)
      deep_arg.xn = xke / Math.pow(a, 1.5)
      /* Long period periodics */ axn = deep_arg.em * Math.cos(deep_arg.omgadf)
      temp = 1 / (a * beta * beta)
      xll = temp * xlcof * axn
      aynl = temp * aycof
      xlt = xl + xll
      ayn = deep_arg.em * Math.sin(deep_arg.omgadf) + aynl
      /* Solve Kepler's Equation */ capu = FMod2p(xlt - deep_arg.xnode)
      temp2 = capu
      i = 0
      var needBreak = false;
      do {
        sinepw = Math.sin(temp2)
        cosepw = Math.cos(temp2)
        temp3 = axn * sinepw
        temp4 = ayn * cosepw
        temp5 = axn * cosepw
        temp6 = ayn * sinepw
        epw = (capu - temp4 + temp3 - temp2) / (1 - temp5 - temp6) + temp2
        if (Math.abs(epw - temp2) <= e6a) needBreak = true //todo: break is not supported
        if (!needBreak) temp2 = epw
      } while ( {
        {
          i += 1;
          i - 1
        } < 10 && !needBreak
      })
      /* Short period preliminary quantities */ ecose = temp5 + temp6
      esine = temp3 - temp4
      elsq = axn * axn + ayn * ayn
      temp = 1 - elsq
      pl = a * temp
      r = a * (1 - ecose)
      temp1 = 1 / r
      rdot = xke * Math.sqrt(a) * esine * temp1
      rfdot = xke * Math.sqrt(pl) * temp1
      temp2 = a * temp1
      betal = Math.sqrt(temp)
      temp3 = 1 / (1 + betal)
      cosu = temp2 * (cosepw - axn + ayn * esine * temp3)
      sinu = temp2 * (sinepw - ayn - axn * esine * temp3)
      u = AcTan(sinu, cosu)
      sin2u = 2 * sinu * cosu
      cos2u = 2 * cosu * cosu - 1
      temp = 1 / pl
      temp1 = ck2 * temp
      temp2 = temp1 * temp
      /* Update for short periodics */ rk = r * (1 - 1.5 * temp2 * betal * x3thm1) + 0.5 * temp1 * x1mth2 * cos2u
      uk = u - 0.25 * temp2 * x7thm1 * sin2u
      xnodek = deep_arg.xnode + 1.5 * temp2 * deep_arg.cosio * sin2u
      xinck = deep_arg.xinc + 1.5 * temp2 * deep_arg.cosio * deep_arg.sinio * cos2u
      rdotk = rdot - deep_arg.xn * temp1 * x1mth2 * sin2u
      rfdotk = rfdot + deep_arg.xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1)
      /* Orientation vectors */ sinuk = Math.sin(uk)
      cosuk = Math.cos(uk)
      sinik = Math.sin(xinck)
      cosik = Math.cos(xinck)
      sinnok = Math.sin(xnodek)
      cosnok = Math.cos(xnodek)
      xmx = -sinnok * cosik
      xmy = cosnok * cosik
      ux = xmx * sinuk + cosnok * cosuk
      uy = xmy * sinuk + sinnok * cosuk
      uz = sinik * sinuk
      vx = xmx * cosuk - cosnok * sinuk
      vy = xmy * cosuk - sinnok * sinuk
      vz = sinik * cosuk
      /* Position and velocity */ pos.x = rk * ux
      pos.y = rk * uy
      pos.z = rk * uz
      vel.x = rdotk * ux + rfdotk * vx
      vel.y = rdotk * uy + rfdotk * vy
      vel.z = rdotk * uz + rfdotk * vz
      /* Phase in radians */ phase = xlt - deep_arg.xnode - deep_arg.omgadf + twopi
      if (phase < 0.0) phase += twopi
      phase = FMod2p(phase)
    }
  }

  def Calculate_User_PosVel(time: Double, geodetic: Predict#geodetic_t, obs_pos: Predict#vector_t, obs_vel: Predict#vector_t): Unit = {
    /* Calculate_User_PosVel() passes the user's geodetic position
        and the time of interest and returns the ECI position and
        velocity of the observer.  The velocity calculation assumes
        the geodetic position is stationary relative to the earth's
        surface. *//* Reference:  The 1992 Astronomical Almanac, page K11. */ var c = .0
    var sq = .0
    var achcp = .0
    geodetic.theta = FMod2p(ThetaG_JD(time) + geodetic.lon)
    /* LMST */ c = 1 / Math.sqrt(1 + f * (f - 2) * Sqr(Math.sin(geodetic.lat)))
    sq = Sqr(1 - f) * c
    achcp = (xkmper * c + geodetic.alt) * Math.cos(geodetic.lat)
    obs_pos.x = achcp * Math.cos(geodetic.theta)
    /* kilometers */ obs_pos.y = achcp * Math.sin(geodetic.theta)
    obs_pos.z = (xkmper * sq + geodetic.alt) * Math.sin(geodetic.lat)
    obs_vel.x = -mfactor * obs_pos.y
    /* kilometers/second */ obs_vel.y = mfactor * obs_pos.x
    obs_vel.z = 0
    Magnitude(obs_pos)
    Magnitude(obs_vel)
  }

  def Calculate_LatLonAlt(time: Double, pos: Predict#vector_t, geodetic: Predict#geodetic_t): Unit = {
    /* Procedure Calculate_LatLonAlt will calculate the geodetic  *//* position of an object given its ECI position pos and time. *//* It is intended to be used to determine the ground track of *//* a satellite.  The calculations  assume the earth to be an  *//* oblate spheroid as defined in WGS '72.                     *//* Reference:  The 1992 Astronomical Almanac, page K12. */ var r = .0
    var e2 = .0
    var phi = .0
    var c = .0
    geodetic.theta = AcTan(pos.y, pos.x)
    /* radians */ geodetic.lon = FMod2p(geodetic.theta - ThetaG_JD(time))
    r = Math.sqrt(Sqr(pos.x) + Sqr(pos.y))
    e2 = f * (2 - f)
    geodetic.lat = AcTan(pos.z, r)
    do {
      phi = geodetic.lat
      c = 1 / Math.sqrt(1 - e2 * Sqr(Math.sin(phi)))
      geodetic.lat = AcTan(pos.z + xkmper * c * e2 * Math.sin(phi), r)
    } while ( {
      Math.abs(geodetic.lat - phi) >= 1E-10
    })
    geodetic.alt = r / Math.cos(geodetic.lat) - xkmper * c
    if (geodetic.lat > pio2) geodetic.lat -= twopi
  }

  def Calculate_Obs(time: Double, pos: vector_t, vel: vector_t, geodetic: geodetic_t, obs_set: vector_t): Unit = {
    /* The procedures Calculate_Obs and Calculate_RADec calculate         *//* the *topocentric* coordinates of the object with ECI position,     *//* {pos}, and velocity, {vel}, from location {geodetic} at {time}.    *//* The {obs_set} returned for Calculate_Obs consists of azimuth,      *//* elevation, range, and range rate (in that order) with units of     *//* radians, radians, kilometers, and kilometers/second, respectively. *//* The WGS '72 geoid is used and the effect of atmospheric refraction *//* (under standard temperature and pressure) is incorporated into the *//* elevation calculation; the effect of atmospheric refraction on     *//* range and range rate has not yet been quantified.                  *//* The {obs_set} for Calculate_RADec consists of right ascension and  *//* declination (in that order) in radians.  Again, calculations are   *//* based on *topocentric* position using the WGS '72 geoid and        *//* incorporating atmospheric refraction.                              */ var sin_lat = .0
    var cos_lat = .0
    var sin_theta = .0
    var cos_theta = .0
    var el = .0
    var azim = .0
    var top_s = .0
    var top_e = .0
    var top_z = .0
    val obs_pos = new vector_t
    val obs_vel = new vector_t
    val range = new vector_t
    val rgvel = new vector_t
    Calculate_User_PosVel(time, geodetic, obs_pos, obs_vel)
    range.x = pos.x - obs_pos.x
    range.y = pos.y - obs_pos.y
    range.z = pos.z - obs_pos.z
    /* Save these values globally for calculating squint angles later... */ rx = range.x
    ry = range.y
    rz = range.z
    rgvel.x = vel.x - obs_vel.x
    rgvel.y = vel.y - obs_vel.y
    rgvel.z = vel.z - obs_vel.z
    Magnitude(range)
    sin_lat = Math.sin(geodetic.lat)
    cos_lat = Math.cos(geodetic.lat)
    sin_theta = Math.sin(geodetic.theta)
    cos_theta = Math.cos(geodetic.theta)
    top_s = sin_lat * cos_theta * range.x + sin_lat * sin_theta * range.y - cos_lat * range.z
    top_e = -sin_theta * range.x + cos_theta * range.y
    top_z = cos_lat * cos_theta * range.x + cos_lat * sin_theta * range.y + sin_lat * range.z
    azim = Math.atan(-top_e / top_s)
    /* Azimuth */ if (top_s > 0.0) azim = azim + Math.PI
    if (azim < 0.0) azim = azim + twopi
    el = ArcSin(top_z / range.w)
    obs_set.x = azim
    /* Azimuth (radians)   */ obs_set.y = el
    /* Elevation (radians) */ obs_set.z = range.w
    /* Range (kilometers)  *//* Range Rate (kilometers/second) */ obs_set.w = Dot(range, rgvel) / range.w
    /* Corrections for atmospheric refraction *//* Reference:  Astronomical Algorithms by Jean Meeus, pp. 101-104    *//* Correction is meaningless when apparent elevation is below horizon */
    /**
      * * The following adjustment for atmospheric refraction is bypassed **
      */
    /* obs_set.y=obs_set.y+Radians((1.02/tan(Radians(Degrees(el)+10.3/(Degrees(el)+5.11))))/60); */ obs_set.y = el

    /**
      * ** End bypass ***
      */
    if (obs_set.y >= 0.0) SetFlag(VISIBLE_FLAG)
    else {
      obs_set.y = el
      /* Reset to true elevation */ ClearFlag(VISIBLE_FLAG)
    }
  }

  def KepCheck(line1: String, line2: String): Boolean = {
    var x = 0
    var sum1 = 0
    var sum2 = 0
    /* Compute checksum for each line */ x = 0
    sum1 = 0
    sum2 = 0
    while ( {
      x <= 67
    }) {
      sum1 += valval(line1.charAt(x).toInt)
      sum2 += valval(line2.charAt(x).toInt)
      x += 1
    }
    /* Perform a "torture test" on the data */ x = (valval(line1.charAt(68).toInt) ^ (sum1 % 10)) | (valval(line2.charAt(68).toInt) ^ (sum2 % 10)) | (line1.charAt(0) ^ '1') | (line1.charAt(1) ^ ' ') | (line1.charAt(7) ^ 'U') | (line1.charAt(8) ^ ' ') | (line1.charAt(17) ^ ' ') | (line1.charAt(23) ^ '.') | (line1.charAt(32) ^ ' ') | (line1.charAt(34) ^ '.') | (line1.charAt(43) ^ ' ') | (line1.charAt(52) ^ ' ') | (line1.charAt(61) ^ ' ') | (line1.charAt(62) ^ '0') | (line1.charAt(63) ^ ' ') | (line2.charAt(0) ^ '2') | (line2.charAt(1) ^ ' ') | (line2.charAt(7) ^ ' ') | (line2.charAt(11) ^ '.') | (line2.charAt(16) ^ ' ') | (line2.charAt(20) ^ '.') | (line2.charAt(25) ^ ' ') | (line2.charAt(33) ^ ' ') | (line2.charAt(37) ^ '.') | (line2.charAt(42) ^ ' ') | (line2.charAt(46) ^ '.') | (line2.charAt(51) ^ ' ') | (line2.charAt(54) ^ '.') | (line1.charAt(2) ^ line2.charAt(2)) | (line1.charAt(3) ^ line2.charAt(3)) | (line1.charAt(4) ^ line2.charAt(4)) | (line1.charAt(5) ^ line2.charAt(5)) | (line1.charAt(6) ^ line2.charAt(6)) | (if (Character.isDigit(line1.charAt(68))) 0
    else 1) | (if (Character.isDigit(line2.charAt(68))) 0
    else 1) | (if (Character.isDigit(line1.charAt(18))) 0
    else 1) | (if (Character.isDigit(line1.charAt(19))) 0
    else 1) | (if (Character.isDigit(line2.charAt(31))) 0
    else 1) | (if (Character.isDigit(line2.charAt(32))) 0
    else 1)
    x == 0
  }

  import java.io.BufferedReader
  import java.io.FileInputStream
  import java.io.IOException
  import java.io.InputStreamReader

  def InternalUpdate(x: Int): Unit = {
    var tempnum = .0
    sat.designator = sat.line1.split("[ \t]+")(2)
    sat.catnum = sat.line1.substring(2, 7).toLong
    sat.year = sat.line1.substring(18, 20).toInt
    sat.refepoch = sat.line1.substring(20, 32).toDouble
    tempnum = 1.0e-5 * sat.line1.substring(44, 50).toDouble
    sat.nddot6 = tempnum / Math.pow(10.0, sat.line1.charAt(51) - '0')
    tempnum = 1.0e-5 * sat.line1.substring(53, 59).toDouble
    sat.bstar = tempnum / Math.pow(10.0, sat.line1.charAt(60) - '0')
    sat.setnum = sat.line1.substring(64, 68).trim.toLong
    sat.incl = sat.line2.substring(8, 16).toDouble
    sat.raan = sat.line2.substring(17, 25).toDouble
    sat.eccn = 1.0e-07 * sat.line2.substring(26, 33).toDouble
    sat.argper = sat.line2.substring(34, 42).toDouble
    sat.meanan = sat.line2.substring(43, 51).toDouble
    sat.meanmo = sat.line2.substring(52, 63).toDouble
    sat.drag = sat.line1.substring(33, 43).toDouble
    sat.orbitnum = sat.line2.substring(63, 68).toLong
  }

  @throws[IOException]
  def ReadDataFiles(file: File): Unit = {
    var name = ""
    var line1 = ""
    var line2 = ""
    val reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)))
    name = reader.readLine
    line1 = reader.readLine
    line2 = reader.readLine
    reader.close()
    if (KepCheck(line1, line2)) {
      sat.name = name
      sat.line1 = line1
      sat.line2 = line2
      InternalUpdate(0)
    }
    else throw new RuntimeException("Invalid TLE File.")
  }


  def DayNum(mPara: Int, d: Int, yPara: Int): Long = {
    var y = yPara
    var m = mPara
    /* This function calculates the day number from m/d/y. */ var dn = 0L
    var mm = .0
    var yy = .0
    if (m < 3) {
      y -= 1
      m += 12
    }
    if (y < 57) y += 100
    yy = y.toDouble
    mm = m.toDouble
    dn = (Math.floor(365.25 * (yy - 80.0)) - Math.floor(19.0 + yy / 100.0) + Math.floor(4.75 + yy / 400.0) - 16.0).toLong
    dn += d + 30 * m + Math.floor(0.6 * mm - 0.3).toLong
    dn
  }

  def PreCalc(): Unit = {
    /* This function copies TLE data from PREDICT's sat structure
        to the SGP4/SDP4's single dimensioned tle structure, and
        prepares the tracking code for the update. */ tle.sat_name = sat.name
    tle.idesg = sat.designator
    tle.catnr = sat.catnum
    tle.epoch = (1000.0 * sat.year.asInstanceOf[Double]) + sat.refepoch
    tle.xndt2o = sat.drag
    tle.xndd6o = sat.nddot6
    tle.bstar = sat.bstar
    tle.xincl = sat.incl
    tle.xnodeo = sat.raan
    tle.eo = sat.eccn
    tle.omegao = sat.argper
    tle.xmo = sat.meanan
    tle.xno = sat.meanmo
    tle.revnum = sat.orbitnum
    calc_squint = 0
    /* Clear all flags */ ClearFlag(ALL_FLAGS)
    /* Select ephemeris type.  This function will set or clear the
         DEEP_SPACE_EPHEM_FLAG depending on the TLE parameters of the
         satellite.  It will also pre-process tle members for the
         ephemeris functions SGP4 or SDP4, so this function must
         be called each time a new tle set is used. */ select_ephemeris(tle)
  }

  def Calc(): Unit = {
    /* This is the stuff we need to do repetitively while tracking. *//* Zero vector for initializations */ val zero_vector = new vector_t
    /* Satellite position and velocity vectors */ val vel = new vector_t
    val pos = new vector_t
    /* Satellite Az, El, Range, Range rate */ val obs_set = new vector_t
    /* Solar ECI position vector  */ val solar_vector = new vector_t
    /* Solar observed azi and ele vector  */ val solar_set = new vector_t
    /* Satellite's predicted geodetic position */ val sat_geodetic = new geodetic_t
    jul_utc = daynum + 2444238.5
    /* Convert satellite's epoch time to Julian  *//* and calculate time since epoch in minutes */ jul_epoch = Julian_Date_of_Epoch(tle.epoch)
    tsince = (jul_utc - jul_epoch) * xmnpda
    age = jul_utc - jul_epoch
    /* Copy the ephemeris type in use to ephem string. */ if (isFlagSet(DEEP_SPACE_EPHEM_FLAG)) ephem = "SDP4"
    else ephem = "SGP4"
    /* Call NORAD routines according to deep-space flag. */ if (isFlagSet(DEEP_SPACE_EPHEM_FLAG)) SDP4(tsince, tle, pos, vel)
    else SGP4(tsince, tle, pos, vel)
    /* Scale position and velocity vectors to km and km/sec */ Convert_Sat_State(pos, vel)
    /* Calculate velocity of satellite */ Magnitude(vel)
    sat_vel = vel.w

    /**
      * All angles in rads. Distance in km. Velocity in km/s *
      */
    /* Calculate satellite Azi, Ele, Range and Range-rate */ Calculate_Obs(jul_utc, pos, vel, obs_geodetic, obs_set)
    /* Calculate satellite Lat North, Lon East and Alt. */ Calculate_LatLonAlt(jul_utc, pos, sat_geodetic)
    /* Calculate squint angle */ if (calc_squint != 0) squint = Math.acos(-((ax * rx + ay * ry + az * rz)) / obs_set.z) / deg2rad
    /* Calculate solar position and satellite eclipse depth. *//* Also set or clear the satellite eclipsed flag accordingly. */ Calculate_Solar_Position(jul_utc, solar_vector)
    Calculate_Obs(jul_utc, solar_vector, zero_vector, obs_geodetic, solar_set)
    val sat_eclipsed_ret = Sat_Eclipsed(pos, solar_vector, eclipse_depth)
    eclipse_depth = sat_eclipsed_ret._2
    if (sat_eclipsed_ret._1) SetFlag(SAT_ECLIPSED_FLAG)
    else ClearFlag(SAT_ECLIPSED_FLAG)
    if (isFlagSet(SAT_ECLIPSED_FLAG)) {
      sat_sun_status = 0
      /* Eclipse */
    }
    else {
      sat_sun_status = 1
      /* In sunlight */
    }
    /* Convert satellite and solar data */ sat_azi = Degrees(obs_set.x)
    sat_ele = Degrees(obs_set.y)
    sat_range = obs_set.z
    sat_range_rate = obs_set.w
    sat_lat = Degrees(sat_geodetic.lat)
    sat_lon = Degrees(sat_geodetic.lon)
    sat_alt = sat_geodetic.alt
    fk = 12756.33 * Math.acos(xkmper / (xkmper + sat_alt))
    fm = fk / 1.609344
    rv = Math.floor((tle.xno * xmnpda / twopi + age * tle.bstar * ae) * age + tle.xmo / twopi).toLong + tle.revnum
    sun_azi = Degrees(solar_set.x)
    sun_ele = Degrees(solar_set.y)
    irk = Math.rint(sat_range).toLong
    isplat = Math.rint(sat_lat).toInt
    isplong = Math.rint(360.0 - sat_lon).toInt
    iaz = Math.rint(sat_azi).toInt
    iel = Math.rint(sat_ele).toInt
    //char logb[50];
    //sprintf (logb, "%d plus %d is %d", a, b, a+b);
    //printf ("[%s] is a string %d chars long\n",buffer,n);
    //logh("in calc");
    ma256 = Math.rint(256.0 * (phase / twopi)).toInt
    if (sat_sun_status != 0) if (sun_ele <= -12.0 && Math.rint(sat_ele) >= 0.0) findsun = '+'
    else findsun = '*'
    else findsun = ' '
  }

  def AosHappens: Boolean = {
    /* This function returns a 1 if the satellite can ever rise above the
           horizon of the ground station. */ var lin = .0
    var sma = .0
    var apogee = .0
    if (sat.meanmo == 0.0) false
    else {
      lin = sat.incl
      if (lin >= 90.0) lin = 180.0 - lin
      sma = 331.25 * Math.exp(Math.log(1440.0 / sat.meanmo) * (2.0 / 3.0))
      apogee = sma * (1.0 + sat.eccn) - xkmper
      if ((Math.acos(xkmper / (apogee + xkmper)) + (lin * deg2rad)) > Math.abs(qth.stnlat * deg2rad)) true
      else false
    }
  }


  def Decayed(time: Double): Boolean = {
    /* This function returns a 1 if it appears that the
        satellite pointed to by 'x' has decayed at the
        time of 'time'.  If 'time' is 0.0, then the
        current date/time is used. */ val satepoch = DayNum(1, 0, sat.year) + sat.refepoch
    satepoch + ((16.666666 - sat.meanmo) / (10.0 * Math.abs(sat.drag))) < time
  }

  def Geostationary: Boolean = {
    /* This function returns a 1 if the satellite
        appears to be in a geostationary orbit */ Math.abs(sat.meanmo - 1.0027) < 0.0002
  }

  def FindAOS: Double = {
    /* This function finds and returns the time of AOS (aostime). */ aostime = 0.0
    if (AosHappens && !Geostationary && !Decayed(daynum)) {
      Calc
      /* Get the satellite in range */ while ( {
        sat_ele < -1.0
      }) {
        daynum -= 0.00035 * (sat_ele * ((sat_alt / 8400.0) + 0.46) - 2.0)
        Calc
      }
      /* Find AOS */ while ( {
        aostime == 0.0
      }) if (Math.abs(sat_ele) < 0.03) aostime = daynum
      else {
        daynum -= sat_ele * Math.sqrt(sat_alt) / 530000.0
        Calc
      }
    }
    aostime
  }

  def FindLOS: Double = {
    lostime = 0.0
    if (Geostationary && AosHappens && !Decayed(daynum)) {
      Calc
      do {
        daynum += sat_ele * Math.sqrt(sat_alt) / 502500.0
        Calc
        if (Math.abs(sat_ele) < 0.03) lostime = daynum
      } while ( {
        lostime == 0.0
      })
    }
    lostime
  }

  def Predict(argv: Array[String]) = {
    val passDetail = new PassDetail
    val startTime = argv(1).toDouble
    obs_geodetic.lat = argv(2).toDouble * deg2rad
    obs_geodetic.lon = argv(3).toDouble * deg2rad
    obs_geodetic.alt = argv(4).toDouble / 1000.0
    obs_geodetic.theta = 0.0
    val div = argv(5).toDouble
    val quit = 0
    var lastel = 0
    var breakout = 0
    PreCalc
    daynum = startTime
    /* Trap geostationary orbits and passes that cannot occur. */
    if (AosHappens && !Geostationary && !Decayed(daynum)) {
      do {
        daynum = FindAOS
        /* Display the pass */ while ( {
          iel >= 0 && quit == 0
        }) {
          //        System.out.println(daynum + "\t" + sat_ele + "\t" + sat_azi + "\t" + sat_range)
          passDetail.appendPosition(daynum, sat_ele, sat_azi, sat_range)
          lastel = iel
          val deltaDaynum = 1.0 / 24.toDouble / 3600 / div
          daynum += deltaDaynum
          Calc
        }
        if (lastel != 0) {
          daynum = FindLOS
          Calc
          if (calc_squint == 0) {
            //sprintf(string,"      %s%4d %4d  %4d  %4d   %4d   %6ld  %4.0f %c\n",Daynum2String(daynum),iel,iaz,ma256,(io_lat=='N'?+1:-1)*isplat,(io_lon=='W'?isplong:360-isplong),irk,squint,findsun);
          }
          else {
            //sprintf(string,"      %s%4d %4d  %4d  %4d   %4d   %6ld  %6ld %c\n",Daynum2String(daynum),iel,iaz,ma256,(io_lat=='N'?+1:-1)*isplat,(io_lon=='W'?isplong:360-isplong),irk,rv,findsun);
          }
        }
        breakout = 1
      } while ( {
        quit == 0 && breakout == 0 && AosHappens && !Decayed(daynum)
      })
      passDetail
    }
    else throw new RuntimeException("Can not predict.")
  }
}

object Main extends App {
  val predict = new Predict
  predict.ReadDataFiles(new File("predict.tle"))
  val passDetail = predict.Predict(Array[String]("predict", "14051.23633", "32.326", "80.026", "5075.0", "1000"))
  passDetail.reveal(10)

//  val a = 1.248923478923472894712893712893712983712983712893791827389123
//  val b = BigDecimal("1.248923478923472894712893712893712983712983712893791827389123")
//  println(f"$a%3.20f")
//  println(b)
}

class PassPosition(val daynum: Double, val ele: Double, val azi: Double, val range: Double) {
  override def equals(obj: scala.Any): Boolean = {
    if (!obj.isInstanceOf[PassPosition]) return false
    val other = obj.asInstanceOf[PassPosition]
    return near(daynum, other.daynum, 0.000001 / 24 / 3600) && near(ele, other.ele, 0.001) && near(azi, other.azi, 0.001) && near(range, other.range, 0.001)
  }

  private def near(a: Double, b: Double, error: Double) = (a >= (b - error)) && (a <= b + error)
}

class PassDetail {
  val passPositions = new ArrayBuffer[PassPosition]()

  def appendPosition(daynum: Double, ele: Double, azi: Double, range: Double) {
    passPositions += new PassPosition(daynum, ele, azi, range)
  }

  def reveal(num: Int) {
    passPositions.slice(0, num).zip(passPositions.slice(1, num + 1)).foreach(z => {
      val v = (z._1.range - z._2.range) / ((z._1.daynum - z._2.daynum)) / 24 / 3600 * 1000
      println(v)
    })
  }

  override def equals(obj: scala.Any): Boolean = {
    if (!obj.isInstanceOf[PassDetail]) return false
    val other = obj.asInstanceOf[PassDetail]
    if (passPositions.size != other.passPositions.size) return false
    passPositions.zip(other.passPositions).forall(z => z._1 == z._2)
  }
}