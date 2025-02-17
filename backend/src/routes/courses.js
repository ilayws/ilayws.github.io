const express = require('express');
const router = express.Router();
const CourseGraph = require('../models/CourseGraph');
const languageService = require('../services/LanguageService');

const courseGraph = new CourseGraph();

// Get supported languages
router.get('/languages', (req, res) => {
    try {
        res.json({
            supported: languageService.getSupportedLanguages(),
            default: languageService.getDefaultLanguage()
        });
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Get all courses
router.get('/', (req, res) => {
    try {
        const { lang } = req.query;
        const courses = courseGraph.getAllCourses();
        const translatedCourses = languageService.translateCourseList(courses, lang);
        res.json(translatedCourses);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Get course details by ID
router.get('/:id', (req, res) => {
    try {
        const { lang } = req.query;
        const courseId = req.params.id;
        const course = courseGraph.getCourseDetails(courseId);
        
        if (!course) {
            return res.status(404).json({ error: 'Course not found' });
        }

        const courseInfo = {
            id: courseId,
            ...course,
            prereqs: courseGraph.getDirectPrerequisites(courseId),
            coreqs: courseGraph.getCorequisites(courseId),
            dependentCourses: courseGraph.getDependentCourses(courseId),
            prereqChainDepth: courseGraph.getPrereqChainDepth(courseId)
        };

        const translatedCourse = languageService.translateCourse(courseInfo, lang);
        res.json(translatedCourse);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Get all prerequisites for a course
router.get('/:id/prerequisites', (req, res) => {
    try {
        const { id } = req.params;
        const { indirect, lang } = req.query;
        const prereqs = courseGraph.getAllPrerequisites(id, indirect !== 'false');
        const translatedPrereqs = languageService.translateCourseList(
            prereqs.map(prereqId => courseGraph.getCourseDetails(prereqId)),
            lang
        );
        res.json(translatedPrereqs);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Get corequisites for a course
router.get('/:id/corequisites', (req, res) => {
    try {
        const { id } = req.params;
        const { lang } = req.query;
        const coreqs = courseGraph.getCorequisites(id);
        const translatedCoreqs = languageService.translateCourseList(
            coreqs.map(coreqId => courseGraph.getCourseDetails(coreqId)),
            lang
        );
        res.json(translatedCoreqs);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Get dependent courses
router.get('/:id/dependents', (req, res) => {
    try {
        const { id } = req.params;
        const { lang } = req.query;
        const dependents = courseGraph.getDependentCourses(id);
        const translatedDependents = languageService.translateCourseList(
            dependents.map(depId => courseGraph.getCourseDetails(depId)),
            lang
        );
        res.json(translatedDependents);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Get available courses based on completed courses
router.post('/available', (req, res) => {
    try {
        const { completedCourses } = req.body;
        const { lang } = req.query;
        
        if (!Array.isArray(completedCourses)) {
            return res.status(400).json({ error: 'completedCourses must be an array' });
        }
        
        const availableCourses = courseGraph.getAvailableCourses(completedCourses);
        const translatedAvailable = languageService.translateCourseList(
            availableCourses.map(courseId => courseGraph.getCourseDetails(courseId)),
            lang
        );
        res.json(translatedAvailable);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Find shortest path between courses
router.get('/path/:fromId/:toId', (req, res) => {
    try {
        const { fromId, toId } = req.params;
        const { lang } = req.query;
        const path = courseGraph.findShortestPath(fromId, toId);
        
        if (!path) {
            return res.status(404).json({ error: 'No path found between courses' });
        }
        
        const translatedPath = languageService.translateCourseList(
            path.map(courseId => courseGraph.getCourseDetails(courseId)),
            lang
        );
        res.json(translatedPath);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Filter courses by criteria
router.get('/filter', (req, res) => {
    try {
        const { semester, difficulty, category, hasPrereqs, hasCoreqs, lang } = req.query;
        const criteria = {
            semester: semester ? parseInt(semester) : undefined,
            difficulty,
            category,
            hasPrereqs: hasPrereqs === 'true',
            hasCoreqs: hasCoreqs === 'true'
        };
        
        const filteredCourses = courseGraph.filterCourses(criteria);
        const translatedFiltered = languageService.translateCourseList(filteredCourses, lang);
        res.json(translatedFiltered);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

module.exports = router; 