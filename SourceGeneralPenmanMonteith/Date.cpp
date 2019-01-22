// $Id: date.C 438 2006-03-06 10:11:20Z bjg $

// Last edited by Trond Reitan, 20/4-2004
// Replaced form() with sprintf in operator<<

#include "Date.h"
#include "DateTime.h"
#include <fstream>
#include <string.h>

// C and Fortran interface:
extern "C" {

    int monthsize(char *datestr)
    {
        char str[20];

        if (strlen(datestr) < 6)
        {
            return 0;
        }

        strncpy(str, datestr, 6);
        strcpy(str + 6, "01");

        Date dt(str);

        if (!dt.legal())
        {
            return 0;
        }

        int size = (int)dt.noDaysInMonth();

        return size;
    }

    int monthsize_(char *datestr)
    {
        return monthsize(datestr);
    }

} // extern "C"

Date NoDate;

/* day_tab includes the days of the months in year/leap year */
static int day_tab[2][13] =
{
    { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 },
    { 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }
};

static int firstdayinmonthtab[2][13] =
{
    { 0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 },
    { 0, 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 }
};

static const char *month_names[12] = { "Jan", "Feb", "Mar", "Apr", "Mai", "Jun",
                                 "Jul", "Aug", "Sep", "Okt", "Nov", "Des"
                               };

static const char *weekday_name[] = { "Monday", "Tuesday", "Wednesday",
                                "Thursday", "Friday", "Saturday", "Sunday"
                              };
static const char *weekday_navn[] = { "Mandag", "Tirsdag", "Onsdag",
                                "Torsdag", "Fredag", "Lørdag", "Søndag"
                              };
static const char *weekday_name_short[] = { "Mon", "Tue", "Wed",
                                      "Thu", "Fri", "Sat", "Sun"
                                    };
static const char *weekday_navn_kort[] = { "Man", "Tir", "Ons",
                                     "Tor", "Fre", "Lør", "Søn"
                                   };

int leap(yearType y)
{
    return  int(y % 4 == 0 && y % 100 != 0 || (y % 400 == 0 && y % 4000 != 0));
}


int legalDate(dayType dd, monthType mm, yearType yyyy)
{
    // illegal month or date < 1
    if (mm < 1 || mm > 12 || dd < 1 || dd > 31 || yyyy < 1753)
    {
        return 0;
    }

    if (dd > day_tab[leap(yyyy)][mm])
    {
        return 0;    // illegal date
    }
    return 1;
}


julType Date::date2jday(void)
{
    return julnum;
}


// The algorithm used in the following 2 functions was described in
// Communications of the ACM, Vol. 6, No. 8 (Aug.63), p.444.

julType Date::date2jday(dayType d, monthType m, yearType y) const
{
    /* Convert a Gregorian Calender date to its corresponding
    Julian day number. Gregorian Calender started on  14,Sep. 1752.
    This function is not valid before that date.
    NB: WE USE 1.JAN.1753 TO DENOTE ILLEGAL DATE.
    FOR EXAMPLE, IF AN ILLEGAL DATE, 35 JAN 1994
    IS GIVEN TO CREATE A DATE VARIABLE, THEN THE
    VALUE STORED IN THERE IS 1.JAN.1753. */

    if (!legalDate(d, m, y))
    {
        d = 1;
        m = 1;
        y = 1753;
    }

    static unsigned long c = 0, ya = 0;
    if (m > 2)
    {
        m -= 3;
    }
    else
    {
        m += 9;
        y--;
    }
    c = y / 100;
    ya = y - 100 * c;
    return ((146097 * c) >> 2) + ((1461 * ya) >> 2) + (153 * m + 2) / 5 + d + 1721119;
} /* date2jday */


void Date::jday2date(dayType &d, monthType &m, yearType &y) const
{
    /* Convert a Julian day number to its corresponding Gregorian
    Calender date. Gregorian Calender started on Sep. 14, 1752.
    This function is not valid before that date. */

    julType j = julnum - 1721119;
    y = (yearType)(((j << 2) - 1) / 146097);
    j = (j << 2) - 1 - 146097 * y;
    d = (dayType)(j >> 2);
    j = ((d << 2) + 3) / 1461;
    d = (dayType)((d << 2) + 3 - 1461 * j);
    d = (d + 4) >> 2;
    m = (5 * d - 3) / 153;
    d = 5 * d - 3 - 153 * m;
    d = (d + 5) / 5;
    y = (yearType)(100 * y + j);

    if (m < 10)
    {
        m += 3;
    }
    else
    {
        m -= 9;
        y++;
    }
} /* jday2date */


Date::Date()
{
    julnum = 2361331;             // Set to 1 JAN 14, which is NoDate
}


Date::Date(const Date &d)
{
    julnum = d.julnum;
}


Date::Date(yearType y, monthType m, dayType d)
{
    julnum = date2jday(d, m, y);
}


Date::Date(const char *date)
{
    //yyyymmdd
    if (strcmp(date, "null") == 0)
    {
        julnum = date2jday(1, 1, 1753);
    }
    else
    {
        char temp[5];
        char temp1[3];
        static dayType dd = 0;
        static monthType mm = 0;
        static yearType yyyy = 0;

        temp[4] = temp1[2] = '\0';
        yyyy = (yearType)atoi(strncpy(temp, date, 4));
        mm = (monthType)atoi(strncpy(temp1, date + 4, 2));
        dd = (dayType)atoi(strncpy(temp1, date + 6, 2));

        julnum = date2jday(dd, mm, yyyy);
    }
}


void skipDelim(istream &is)
{
    // Help function to read and discard the delimiting character
    char c;
    if (!is.good())
    {
        return;
    }
    is >> c;
    while (is.good() && !isalnum(c))
    {
        is >> c;
    }
    if (is.good())
    {
        is.putback(c);
    }
}


julType Date::parseDate(istream& is) const
{
    /* parse dates of the format yyyy.mm.dd
    or dd.mm.yyyy. Delimiter need not be '.',
    it can be any of the alphanumeric character */
    static dayType d = 0;
    static monthType m = 0;
    static yearType y = 0;
    static unsigned temp = 0;

    if (is.good())
    {
        skipDelim(is);
        is >> temp;
        skipDelim(is);
        if (is.eof())
        {
            return 0;
        }
        if (temp / 100 > 0)
        {
            // date formate is: yyyy.mm.dd
            y = temp;
            is >> m;
            skipDelim(is);
            if (is.eof())
            {
                return 0;
            }
            is >> d;
        }
        else
        {
            // date formate is: dd.mm.yyyy
            d = temp;
            is >> m;
            skipDelim(is);
            if (is.eof())
            {
                return 0;
            }
            is >> y;
        }
    }
    if (!is.good())
    {
        return 0;
    }
    return Date((yearType)y, m, d).julnum;
} // end parseDate


Date::Date(istream &is)
{
    julnum = parseDate(is);
}


ostream& operator<<(ostream& os, Date &d)
{
    if (d == NoDate)
    {
        os << "NULL";
    }
    else
    {
        static dayType dd = 0;
        static monthType mm = 0;
        static yearType yyyy = 0;
        static char buffer[15];

        d.jday2date(dd, mm, yyyy);
        sprintf(buffer, "%02d/%02d/%d ", dd, mm, yyyy);
        os << buffer;
    }

    return os;
}


int Date::Get_weekday_index(void)
{
    static Date monday(1999, 1, 25);

    int weekday, diff = (int)abs((long)(julnum - monday.julnum));
    if (julnum < monday.julnum)
    {
        weekday = (7 - diff % 7) % 7;
    }
    else
    {
        weekday = diff % 7;
    }
    return weekday;
}


WEEKDAY Date::Get_weekday(void)
{
    return (WEEKDAY)Get_weekday_index();
}


const char *Date::Get_weekday_name_eng(void)
{
    return weekday_name[Get_weekday_index()];
}


const char *Date::Get_weekday_name_norw(void)
{
    return weekday_navn[Get_weekday_index()];
}


const char *Date::Get_weekday_name_eng_short(void)
{
    return weekday_name_short[Get_weekday_index()];
}


const char *Date::Get_weekday_name_norw_short(void)
{
    return weekday_navn_kort[Get_weekday_index()];
}


void Date::now()
{
    time_t clk = time(0);
    tm* now;
    now = localtime(&clk);
    if (now->tm_year < 50)
    {
        julnum = date2jday(now->tm_mday, now->tm_mon + 1, now->tm_year + 2000);
    }
    else
    {
        julnum = date2jday(now->tm_mday, now->tm_mon + 1, now->tm_year + 1900);
    }
}


int Date::legal() const
{
    return (*this != NoDate) ? 1 : 0;
}


Date& Date::operator=(const DateTime &d1)
{
    julnum = d1.julnum;
    return *this;
}


yearType Date::getYear() const
{
    static dayType dd = 0;
    static monthType mm = 0;
    static yearType yyyy = 0;
    jday2date(dd, mm, yyyy);
    return yyyy;
}


monthType Date::getMonth() const
{
    static dayType dd = 0;
    static monthType mm = 0;
    static yearType yyyy = 0;
    jday2date(dd, mm, yyyy);
    return mm;
}


dayType Date::getDay() const
{
    static dayType dd = 0;
    static monthType mm = 0;
    static yearType yyyy = 0;
    jday2date(dd, mm, yyyy);
    return dd;
}

int Date::leap() const
{
    return (::leap(getYear())) ? 1 : 0;
}



const char *Date::getMonthString() const
{
    if (*this == NoDate)
    {
        return "null";
    }
    else
    {
        static dayType d = 0;
        static monthType m = 0;
        static yearType y = 0;
        jday2date(d, m, y);
        return month_names[m - 1];
    }
}


char *Date::getChDate() const
{
    static char date[12];

    for (register int i = 0; i < 12; i++)
    {
        date[i] = 0;
    }
    if (*this == NoDate)
    {
        sprintf(date, "null");
    }
    else
    {
        static dayType d = 0;
        static monthType m = 0;
        static yearType y = 0;

        jday2date(d, m, y);

        sprintf(date, "%d %s %02d", y, getMonthString(), d);
    }

    return date;
}


int Date::dayOfYear() const
{
    static dayType dd = 0;
    static monthType mm = 0;
    static yearType yy = 0;
    jday2date(dd, mm, yy);

    return firstdayinmonthtab[::leap(yy)][mm] + dd;
}


int Date::dayIndex() const
{
    static dayType dd = 0;
    static monthType mm = 0;
    static yearType yy = 0;
    jday2date(dd, mm, yy);

    if (mm == 2 && dd == 29)
    {
        return -1;
    }

    return firstdayinmonthtab[0][mm] + dd;
}

int Date::daysInYear() const
{
    return  (leap()) ? 366 : 365;
}


int Date::firstDayOfMonth() const
{
    static dayType dd = 0;
    static monthType mm = 0;
    static yearType yyyy = 0;
    register int d = 0, i;
    jday2date(dd, mm, yyyy);
    for (i = 1; i < mm; i++)
    {
        d += day_tab[leap()][i];
    }
    d++;
    return d;
}


int Date::noDaysInMonth() const
{
    static dayType dd = 0;
    static monthType mm = 0;
    static yearType yyyy = 0;
    jday2date(dd, mm, yyyy);
    return day_tab[leap()][mm];
}



char *Date::charDate() const
{
    char *pt = new char[10];
    if (*this == NoDate)
    {
        sprintf(pt, "null");
    }
    else
    {
        sprintf(pt, "%d/%d/%d", getYear(), getDay(), getMonth());
    }
    return (pt);
}


int Date::idxOfYear() const
{
    static dayType dd = 0;
    static monthType mm = 0;
    static yearType yy = 0;
    jday2date(dd, mm, yy);

    return firstdayinmonthtab[::leap(yy)][mm] + dd - 1;
}


int Date::getDayOfWeek() const
{
    // 15. sept. 1752 was a Monday (day 1), so:
    return (julnum - date2jday(15, 9, 1752) + 1) % 7;
}


int overlap_date(const Date &start1, const Date &end1,
                 const Date &start2, const Date &end2)
{
    if (start1 != NoDate && end2 != NoDate && start1 > end2)
    {
        return 0;
    }

    if (end1 != NoDate && start2 != NoDate && end1 < start2)
    {
        return 0;
    }

    return 1;
}

#ifdef MAIN

int main(int argc, char **argv)
{
    //cout << "hei, hei!" << endl;

    Date d((yearType)1, 1, 1752), dd;
    cout << d << endl;
    cout << d.date2jday(1, 1, 1753) << endl;
    cout << dd << endl;
    Date d2(1998, 1, 12);
    cout << "12/1-1998:" << d2.Get_weekday_name_norw() << endl;
    cout << d2.idxOfYear() << endl;

    return 0;
}

#endif
