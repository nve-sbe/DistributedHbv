// $Id: date_time.h 442 2006-03-07 07:58:51Z bjg $
#ifndef date_time_h
#define date_time_h
//
//      Class name : DateTime
//
//      Description : This is a subclass of Date and Time, which represents
//        a date-time by a julian day for date and total seconds for time
//        which are private.
//
// Altered by Erik Taraldsen Thu Jul  1 11:14:23 CEST 1999

// ********************************************
//                   License
// This file and it's associated C-file is
// licensed as LGPL. For a description of
// this, see the license.txt file or
// http://www.gnu.org/licenses/licenses.html
// ********************************************

#include "Date.h"
#include "CTime.h"
#include <string.h>
#include <stdio.h> //erikt Added to avoid implisit declaration


//extern int sscanf( const char *str, const char *format, ...); //erikt

int legalDateTime(yearType y, monthType m, dayType d, short h, short min,
                  short s);

class DateTime : public  Date, public CTime
{


    DateTime(julType j, totalSecondsType ts) : Date(j), CTime(ts) { };

public:
    DateTime();
    DateTime(const char *date, const char* time);
    DateTime(yearType yy, monthType mm = 1, dayType dd = 1, hourType hh = 0,
             minuteType mn = 0, secondType second = 0);

    DateTime(const Date &d);
    DateTime(const DateTime& d);
    DateTime(const char *datetime);

    DateTime(char *datetime, int form); // se syCh(int form)

    void now();                         // sets the datetime to now.
    int legal() const;          // checks for legal datetime.

    DateTime Floor(int precision_min); // round down
    // to the given precision (in minutes)
    DateTime Ceiling(int precision_min); // round up
    // to the given precision (in minutes)

    // difference in minutes between this datetime and the second,
    // given datetime
    long int difference_minutes(DateTime &dt2);

    DateTime &add_minutes(long int num_minutes);

    // assignment operator
    DateTime& operator=(const DateTime &dt2);
    DateTime& operator=(const Date &dt2);

    // logical operators
    friend int operator== (const DateTime &dt1, const DateTime &dt2);
    friend int operator!= (const DateTime &dt1, const DateTime &dt2);
    friend int operator<  (const DateTime &dt1, const DateTime &dt2);
    friend int operator<= (const DateTime &dt1, const DateTime &dt2);
    friend int operator>  (const DateTime &dt1, const DateTime &dt2);
    friend int operator>= (const DateTime &dt1, const DateTime &dt2);

    // arithmetical operators

    friend DateTime operator+ (const DateTime &dt, int sec);
    friend DateTime operator+ (int sec, const DateTime &dt);

    friend unsigned long operator- (const DateTime &dt1, const DateTime &dt2);
    // this func. returns the absolute difference between the two times
    // in seconds

    friend DateTime operator- (const DateTime &d, int sec);
    friend DateTime operator- (int sec, const DateTime &d);
    // this functions substract 'sec' seconds from the datetime 'd'
    // and returns the substracted datetime

    friend void operator-= (DateTime &dt, int sec);
    friend void operator+= (DateTime &dt, int sec);

    // others

    // d1 <= dt2 && dt2 <= dt3
    friend int between(const DateTime &d1, const DateTime &dt2,
                       const DateTime &dt3);

    // Does two periods overlap?
    friend int overlap(const DateTime &start1, const DateTime &end1,
                       const DateTime &start2, const DateTime &end2);

    friend ostream& operator<<(ostream&, const DateTime &);
    friend istream& operator>>(istream&, DateTime &);
    void Print();
    void getDateTime(yearType &yyyy, monthType &mm, dayType &dd,
                     hourType &hh, minuteType &mn, secondType &ss);

    char *getChDateTime() const;

    char *getRfc822Str(int UT_offset = 1) const;

    char *syCh(int form = 0) const;             // Sybase format
    // 0 - mm/dd/yyyy hh:mm:ss
    // 1 - dd/mm/yyyy hh:mm
    // 2 - yyyymmddhhmm
    // 3 - yyyy-mm/dd hh:mm
    // 4 - dd/mm/yyyy



    void internal(char str[13]);

    int  isNull() const;


    DateTime StartOfDay() const;                // yyyymmdd 00:00:00
    DateTime EndOfDay() const;                  // yyyymmdd 23:59:00

    DateTime StartOfMonth() const;              // yyyymmdd 00:00:00
    DateTime EndOfMonth() const;                // yyyymmdd 23:59:0

    DateTime StartOfYear() const;               // yyyymmdd 00:00:00
    DateTime EndOfYear() const;                 // yyyymmdd 23:59:00
    double as_floating_point_year(void);


    // Conversion
    //operator Date() { return (Date&)(*this);};
};

//#ifndef LINUX  // erikt If we are not on linux NoDateTime is defined in this way
//extern DateTime NoDateTime;
//#endif
//#ifdef LINUX  // erikt


static DateTime NoDateTime(1753, 1, 1, 0, 0, 0); //erikt
//#endif

inline int operator== (const DateTime &dt1, const DateTime &dt2)
{
    return int(dt1.julnum == dt2.julnum && dt1.tsec == dt2.tsec);
}

inline int operator!= (const DateTime &dt1, const DateTime &dt2)
{
    return int(dt1.julnum != dt2.julnum || dt1.tsec != dt2.tsec);
}


inline int operator> (const DateTime &dt1, const DateTime &dt2)
{
    return dt2 < dt1;
}

inline int operator>= (const DateTime &dt1, const DateTime &dt2)
{
    return dt2 <= dt1;
}

inline DateTime maxdt(DateTime dt1, DateTime dt2)
{
    if (dt1 != NoDateTime && (dt1 > dt2 || dt2 == NoDateTime))
    {
        return dt1;
    }
    else
    {
        return dt2;
    }
}

inline DateTime mindt(DateTime dt1, DateTime dt2)
{
    if (dt1 != NoDateTime && (dt1 < dt2 || dt2 == NoDateTime))
    {
        return dt1;
    }
    else
    {
        return dt2;
    }
}

inline int between(const DateTime &dt1, const DateTime &dt2,
                   const DateTime &dt3)
{
    return int((dt2 != NoDateTime) && ((Date)dt2 != NoDate) &&
               (dt1 == NoDateTime || dt1 <= dt2) &&
               (dt3 == NoDateTime || dt2 <= dt3));
}


extern const char *Month[12];
extern const char *Month_eng[12];
extern const char *Month_short[12];
extern const char *Month_eng_short[12];


// take care of converting operations regarding the
// special DateTime format
class Datetimes
{
public:
    char timestr[20];

    const char* Dt2Str(DateTime dt, int type); // returns a converted datetime
    DateTime Str2Dt(char *dtstr); // returns a converted string
};


#endif
