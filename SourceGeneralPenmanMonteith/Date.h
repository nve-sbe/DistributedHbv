// $Id: date.h 320 2006-01-17 13:27:29Z bjg $
//      File name   : date.h
//      Class name  : Date
//
//      Description : Describes a date and the operations belongs to it.
//      A date is represented as a julian day number (private member of
//      the class). Earliest date that can be represented is 14 Sep. 1752,
//      which we use to denote illegal date.
//
// Alterd by Erik Taraldsen Thu Jul  1 10:54:33 CEST 1999

// ********************************************
//                   License
// This file and it's associated C-file is
// licensed as LGPL. For a description of
// this, see the license.txt file or
// http://www.gnu.org/licenses/licenses.html
// ********************************************

#ifndef date_h
#define date_h

#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <stdio.h> //erikt Added to avoid implisit declaration
#include <iostream>
#include <fstream>

using namespace std;

class DateTime;

typedef unsigned short dayType;
typedef unsigned short monthType;
typedef unsigned short yearType;
typedef unsigned long julType;
enum WEEKDAY
{
    MONDAY = 0, TUESDAY = 1, WEDNESDAY = 2, THURSDAY = 3,
    FRIDAY = 4, SATURDAY = 5, SUNDAY = 6
};

int leap(yearType y);
int legalDate(dayType dd, monthType mm, yearType yyyy);


// extern int sprintf(char *str, const char *format, ...);

class Date
{
    friend class DateTime;
protected:
    julType julnum;
    Date(julType j)
    {
        julnum = j;
    };

public:
    julType date2jday(void);
    julType date2jday(dayType, monthType, yearType) const;
    /* Convert a Gregorian Calender date to its corresponding
    Julian day number. Gregorian Calender started on Sep. 14, 1752
    (we use this to denote an illegal date).
    This function is not valid before that date. */

    void jday2date(dayType&, monthType&, yearType&) const;
    /* Convert a Julian day number to its corresponding Gregorian
    Calender date. Gregorian Calender started on Sep. 14, 1752.
    This function is not valid before that date. */

public:
    Date();
    Date(istream& is);
    Date(const char *);
    Date(const Date &);
    Date(yearType y, monthType m, dayType d = 1);

    void now(); // sets the date to today
    int legal() const; /* returns true if a legal date.
					   OBS: illegal date is stored as 1752.sep.14*/

    // assignment operator
    Date& operator=(const Date &d1)
    {
        julnum = d1.julnum;
        return *this;
    };
    Date& operator=(const DateTime &d1);

    // logical operators
    friend int operator== (const Date &d1, const Date &d2)
    {
        return int(d1.julnum == d2.julnum);
    };
    friend int operator!= (const Date &d1, const Date &d2)
    {
        return int(d1.julnum != d2.julnum);
    };
    friend int operator<(const Date &d1, const Date &d2)
    {
        return int(d1.julnum < d2.julnum);
    };
    friend int operator<=(const Date &d1, const Date &d2)
    {
        return int(d1.julnum <= d2.julnum);
    };
    friend int operator>(const Date &d1, const Date &d2)
    {
        return int(d2 < d1);
    };
    friend int operator>= (const Date &d1, const Date &d2)
    {
        return int(d2 <= d1);
    };

    // arithmetical operators

    // adds some days 'dd' to a date and returns the resulting date
    friend Date operator+ (const Date &d, int dd)
    {
        return Date(d.julnum + dd);
    };
    friend Date operator+ (int dd, const Date &d)
    {
        return Date(dd + d.julnum);
    };

    //adds days 'dd' to the given date 'd';
    //note: if you don't want the original date changed, use + mentioned above
    friend void operator+= (Date &d, int dd)
    {
        d.julnum += dd;
    };

    //returns the difference between 2 dates in number of days
    friend int operator- (const Date &d1, const Date &d2)
    {
        return (d1.julnum - d2.julnum);
    };

    // substract 'dd' days from a date 'd' and returns the resulting date
    friend Date operator- (const Date &d, int dd)
    {
        return Date(d.julnum - dd);
    };

    friend void operator-= (Date &d, int dd)
    {
        d.julnum -= dd;
    };

    // others
    friend int between(const Date &d3, const Date &d1, const Date &d2)
    {
        return int(d3.julnum >= d1.julnum && d3.julnum <= d2.julnum);
    };

    friend int overlap_date(const Date &start1, const Date &end1,
                            const Date &start2, const Date &end2);

    friend ostream& operator<<(ostream&, Date &);
    //  void print();

    char *getChDate() const;              // yyyy mmm dd
    char *charDate() const;  // yyyy:mm:dd
    yearType getYear() const;
    monthType getMonth() const;
    dayType getDay() const;
    const char *getMonthString() const;
    int getDayOfWeek() const; // sunday = 0, saturday = 6

    int leap() const; // whether the current date-object's year is leap

    julType parseDate(istream&) const; // help function for Date(istream&)

    int dayOfYear() const;    // returns the corresponding day of the year
    // starting at 1


    int dayIndex() const;   // returns the corresponding day of the year
    // starting at 1, ingoring leap days


    int daysInYear() const; // returns no. of days in the year: 365 or 366
    int firstDayOfMonth() const;  // returns the first day of the month
    int noDaysInMonth() const;    // returns no. of days in the month

    int idxOfYear() const;                // return pos in year starting at

    int Get_weekday_index(void);
    WEEKDAY Get_weekday(void);
    const char *Get_weekday_name_eng(void);
    const char *Get_weekday_name_norw(void);
    const char *Get_weekday_name_eng_short(void);
    const char *Get_weekday_name_norw_short(void);
};

extern Date NoDate;

#endif
