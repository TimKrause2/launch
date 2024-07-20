#ifndef EVENT_H
#define EVENT_H


#include <list>
#include <iostream>
#include <stdio.h>

struct TimeSpec
{
    char str_buff[64];

    long sec;
    long nsec;
    TimeSpec()
    {
        sec = 0;
        nsec = 0;
    }
    TimeSpec(long sec, long nsec):
        sec(sec),
        nsec(nsec)
    {
    }
    TimeSpec(double t)
    {
        sec = (long)t;
        nsec = (long)((t-sec)*1e9);
    }
    void print(void)
    {
        std::cout << "sec:" << sec << " nsec:" << nsec << std::endl;
    }
    char *to_str(void)
    {
        sprintf(str_buff,"%ld.%03ld",sec,labs(nsec/1000000));
        return str_buff;
    }
    double to_double(void)
    {
        return (double)sec + (double)nsec*1e-9;
    }
    TimeSpec operator=(TimeSpec rhs)
    {
        sec = rhs.sec;
        nsec = rhs.nsec;
        return *this;
    }
    TimeSpec operator+(TimeSpec ts)
    {
        long new_sec = sec + ts.sec;
        long new_nsec = nsec + ts.nsec;
        // check for nsec overflow
        if(new_nsec >= 1000000000){
            new_sec++;
            new_nsec -= 1000000000;
        }else if(new_nsec <=-1000000000){
            new_sec--;
            new_nsec += 1000000000;
        }
        // check for sign mismatch
        if(new_sec > 0){
            if(new_nsec < 0){
                new_sec--;
                new_nsec += 1000000000;
            }
        }else if(new_sec < 0){
            if(new_nsec > 0){
                new_sec++;
                new_nsec -= 1000000000;
            }
        }
        return TimeSpec(new_sec, new_nsec);
    }
    TimeSpec operator+=(TimeSpec ts)
    {
        *this = *this + ts;
        return *this;
    }
    TimeSpec operator+=(double dt)
    {
        *this = *this + TimeSpec(dt);
        return *this;
    }
    TimeSpec operator+(double t)
    {
        TimeSpec rhs(t);
        return *this + rhs;
    }
    TimeSpec operator-(TimeSpec ts)
    {
        TimeSpec rhs(-ts.sec, -ts.nsec);
        return *this + rhs;
    }
    TimeSpec operator-(double t)
    {
        TimeSpec rhs(-t);
        return *this + rhs;
    }
    bool operator>(TimeSpec ts)
    {
        if(sec>ts.sec){
            return true;
        }else if(sec==ts.sec){
            if(nsec>ts.nsec){
                return true;
            }
        }
        return false;
    }
    bool operator>=(TimeSpec ts)
    {
        if(sec>ts.sec){
            return true;
        }else if(sec==ts.sec){
            if(nsec>=ts.nsec){
                return true;
            }
        }
        return false;
    }
    bool operator<(TimeSpec ts)
    {
        if(sec<ts.sec){
            return true;
        }else if(sec==ts.sec){
            if(nsec<ts.nsec){
                return true;
            }
        }
        return false;
    }
    bool operator<=(TimeSpec ts)
    {
        if(sec<ts.sec){
            return true;
        }else if(sec == ts.sec){
            if(nsec <= ts.nsec){
                return true;
            }
        }
        return false;
    }
};

class Event
{
public:
    TimeSpec t;
    Event(TimeSpec t_absolute);
    virtual ~Event();
    bool Expired(TimeSpec &tnow);
    virtual void Action(void)=0;
};

class Satellite;

class SatelliteEvent : public Event
{
    Satellite *satellite;
public:
    SatelliteEvent(Satellite* satellite, TimeSpec t_absolute);
    ~SatelliteEvent();
    void Action(void);
};

class IntegratorEvent : public Event
{
public:
    IntegratorEvent(TimeSpec t_absolute);
    void Action(void){}
};

class EventQueue
{
    TimeSpec tnow;
    std::list<Event*> events;
public:
    void AddEvent(Event* event);
    bool Empty(void);
    void Check(void);
    void Update(double dt);
    double NextEvent(void);
    TimeSpec AbsoluteTime(double time_from_now);
};

#endif
