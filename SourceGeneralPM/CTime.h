// $Id: ctime.h 438 2006-03-06 10:11:20Z bjg $
//      File: ctime.h
//      Class name : Time
//
//      Description : Describes a time class and the operators
//      Time is represented as total seconds (protected member)
//      Max time that can be stored is 2**sizeof(int)-1.

// ********************************************
//                   License
// This file and it's associated C-file is
// licensed as LGPL. For a description of
// this, see the license.txt file or
// http://www.gnu.org/licenses/licenses.html
// ********************************************

#ifndef ctime_h
#define ctime_h

#include <iostream>
#include <fstream>
#include <time.h>

using namespace std;

typedef unsigned short hourType;
typedef unsigned short minuteType;
typedef unsigned short secondType;
typedef unsigned int totalSecondsType;

int legalTime(const short h, const short m, const short s);

class CTime
{

protected:
    totalSecondsType tsec;
    CTime(totalSecondsType t) : tsec(t) {};

public:
    void totSec2time(hourType&, minuteType&, secondType&) const;

    /* Convert seconds to hhmmss CTime (hour:minute:sec). This function is
    not valid for CTime whose total seconds is greter
    than 2**sizeof(int)-1. */

    totalSecondsType time2totSec(hourType, minuteType, secondType) const;

    /* Convert an hhmmss CTime (hour:minute:sec) to seconds.
    This function is not valid for CTime whose total seconds is greter
    than 2**sizeof(int)-1. */

public:
    CTime() {};
    CTime(istream&);
    CTime(const CTime &);
    CTime(const char *);
    CTime(hourType, minuteType, secondType = 0);

    void now();
    int legal() const;

    // assignment operator
    CTime& operator= (const CTime &t2)
    {
        tsec = t2.tsec;
        return *this;
    };

    // assignment operators
    friend int operator== (const CTime &t1, const CTime &t2)
    {
        return int(t1.tsec == t2.tsec);
    };
    friend int operator!= (const CTime &t1, const CTime &t2)
    {
        return int(t1.tsec != t2.tsec);
    };
    friend int operator<(const CTime &t1, const CTime &t2)
    {
        return int(t1.tsec < t2.tsec);
    };
    friend int operator<=(const CTime &t1, const CTime &t2)
    {
        return int(t1.tsec <= t2.tsec);
    };
    friend int operator>(const CTime &t1, const CTime &t2)
    {
        return t2 < t1;
    };
    friend int operator>= (const CTime &t1, const CTime &t2)
    {
        return t2 <= t1;
    };

    friend CTime operator+ (const CTime &t, int sec)
    {
        return CTime(t.tsec + sec);
    };
    friend CTime operator+ (const CTime &t1, const CTime &t2)
    {
        return CTime(t1.tsec + t2.tsec);
    };
    friend void operator+= (CTime &t, int sec)
    {
        t.tsec += sec;
    };

    friend int operator- (const CTime &t1, const CTime &t2)
    {
        return (int)(t1.tsec - t2.tsec);
    };
    friend CTime operator- (CTime &t, int sec)
    {
        return CTime(t.tsec - sec);
    };
    friend void operator-= (CTime &t, int sec)
    {
        t.tsec = t.tsec - sec;
    };
    // friend void operator-= (CTime &t, const CTime &t2)
    //    { t.tsec = t.tsec -t2.tsec; };

    friend int between(const CTime &t3, const CTime &t1, const CTime &t2)
    {
        return int(t3.tsec >= t1.tsec && t3.tsec <= t2.tsec);
    };

    friend ostream& operator<<(ostream&, const CTime &);
    //void print();

    hourType getHour() const;
    minuteType getMinute() const;
    secondType getSecond() const;
    char *getChTime() const;

    int minofday() const;               // [0 - 1399]

    // help function for CTime(istream&)
    totalSecondsType parseTime(istream&)const;
};
#endif
